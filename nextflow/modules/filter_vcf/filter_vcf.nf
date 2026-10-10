nextflow.enable.types = true

process FILTER_VCF {

      tag "${interval_id}"

      conda "${moduleDir}/environment.yml"

      input:
      record(
            interval_id: String,
            interval_bed: Path,
            interval_index: Path,
            variants: Path,
            variants_index: Path
      )
      mask_bed: Path
      popmap: Path
      missing_summary: Path

      output:
      record(
            interval_id: interval_id,
            interval_bed: interval_bed,
            interval_index: interval_index,

            all_bcf: file("${interval_id}.filt.bcf"),
            all_bcf_index: file("${interval_id}.filt.bcf.csi"),

            snp_bcf: file("${interval_id}.snp.bcf"),
            snp_bcf_index: file("${interval_id}.snp.bcf.csi"),

            indel_bcf: file("${interval_id}.indel.bcf"),
            indel_bcf_index: file("${interval_id}.indel.bcf.csi"),

            invariant_bcf: file("${interval_id}.invariant.bcf"),
            invariant_bcf_index: file("${interval_id}.invariant.bcf.csi"),

            filter_hist: file("${interval_id}.filter_hist.tsv"),
            mask_hist: file("${interval_id}.mask_hist.tsv"),
            samples_to_keep: file("${interval_id}.samples.txt")
      )

      script:
      """
      #!/usr/bin/env bash
      set -euo pipefail
      
      # Source dependent functions
      source "\$(command -v functions.sh)"

      ###########################################
      # Set up mask filters
      ###########################################

      case "${mask_bed}" in
            *.gz) gzip -dc "${mask_bed}" ;;
            *)    cat "${mask_bed}" ;;
      esac |
      awk '
            BEGIN { FS=OFS="\\t" }
            /^#/ || /^track[[:space:]]/ || /^browser[[:space:]]/ || NF == 0 { next }
            {
                  sub(/\\r\$/, "", \$NF)

                  if (NF < 3 || \$2 !~ /^[0-9]+\$/ ||
                      \$3 !~ /^[0-9]+\$/ || (\$3+0) <= (\$2+0)) {
                        print "ERROR: invalid mask BED line " NR > "/dev/stderr"
                        exit 1
                  }

                  name = (NF >= 4 && \$4 != "" && \$4 != "." ? \$4 : "Excluded")
                  if (name != "Included")
                        print \$1, \$2, \$3, name
            }
      ' |
      LC_ALL=C sort -u -k1,1 -k2,2n -k3,3n -k4,4 > vcf_masks.bed

      mkdir -p mask_beds
      cut -f4 vcf_masks.bed | LC_ALL=C sort -u > mask_names.txt

      printf 'filter_id\\tmask_name\\n' > "${interval_id}.mask_labels.tsv"
      MASK_FILTERS=()

      while IFS= read -r mask_name; do
            slug=\$(printf '%s' "\$mask_name" | LC_ALL=C tr -c 'A-Za-z0-9_' '_')
            filter_id="MASK_\${slug}_FAIL"
            bed_file="mask_beds/\${filter_id}.bed"

            if [[ -e "\$bed_file" ]]; then
                  echo "ERROR: mask names produce the same label: \$filter_id" >&2
                  exit 1
            fi

            MASK_NAME="\$mask_name" awk '
                  BEGIN { FS=OFS="\\t" }
                  \$4 == ENVIRON["MASK_NAME"] { print \$1, \$2, \$3 }
            ' vcf_masks.bed > "\$bed_file"

            MASK_FILTERS+=("\$filter_id")
            printf '%s\\t%s\\n' "\$filter_id" "\$mask_name" \
                  >> "${interval_id}.mask_labels.tsv"
      done < mask_names.txt

      apply_mask_filters() {
            local index="\${1:-0}"

            if (( index == \${#MASK_FILTERS[@]} )); then
                  cat
            else
                  bcftools filter -Ou \
                        -s "\${MASK_FILTERS[\$index]}" \
                        -m+ --mask-overlap 1 \
                        -M "mask_beds/\${MASK_FILTERS[\$index]}.bed" |
                  apply_mask_filters "\$((index + 1))"
            fi
      }

      ###########################################
      # Set up per-sample missing filters
      ###########################################

      # Find samples above the missing fraction filter
      awk -v thr="${params.sample_max_missing}" 'NR==1 {next} \$4!="NA" && (\$4+0) < thr {print \$1}' "${missing_summary}" > ${interval_id}.samples.txt

      ###########################################
      # Set up Population filters
      ###########################################

      # Create sample_groups.tsv:
      # first keep only samples in samples.txt
      # then drop populations with fewer than MIN_SAMPLES_PER_POP retained samples
      awk -v n="${params.vcf_population_min_samples}" '
      BEGIN { FS=OFS="\\t" }
      NR==FNR {
            keep[\$1] = 1
            next
      }
      (\$1 in keep) {
            pop_name = \$2
            gsub(/[[:space:]_]+/, "", pop_name)

            count[pop_name]++
            sample[NR] = \$1
            pop[NR]    = pop_name
      }
      END {
            for (i = 1; i <= NR; i++) {
                  if (sample[i] != "" && count[pop[i]] >= n) {
                  print sample[i], pop[i]
                  }
            }
      } 
      ' "${interval_id}.samples.txt" "${popmap}" > sample_groups.tsv

      # Return success only for enabled threshold values.
      enabled() {
            local value="\${1:-}"

            case "\${value,,}" in
                  ""|null|na)
                  return 1
                  ;;
                  *)
                  return 0
                  ;;
            esac
      }

      join_by() {
            local separator="\$1"
            shift

            local value
            local output=""

            for value in "\$@"; do
                  [[ -n "\$value" ]] || continue
                  [[ -z "\$output" ]] || output+="\$separator"
                  output+="\$value"
            done

            printf '%s' "\$output"
      }

      case "${params.vcf_population_fail_mode}" in
            ALL)
                  population_join=" && "
                  ;;
            ANY)
                  population_join=" || "
                  ;;
            *)
                  echo "ERROR: vcf_population_fail_mode must be ALL or ANY, got '${params.vcf_population_fail_mode}'" >&2
                  exit 1
                  ;;
      esac

      mapfile -t POPS < <(
            cut -f2 "sample_groups.tsv" |
                  tr ',' '\\n' |
                  awk 'NF' |
                  LC_ALL=C sort -u
      )

      if (( \${#POPS[@]} == 0 )); then
            echo "ERROR: no populations found in sample_groups.tsv" >&2
            exit 1
      fi

      # Build either a global or per-population bcftools exclusion expression.
      #
      # Examples:
      #
      #   make_filter_expr DP "<" GLOBAL "SNP:5" "INDEL:10"
      #   make_filter_expr MAF "<" POP "SNP:0.05" "INDEL:0.05"
      #
      # Disabled thresholds may be blank, null or NA.
      make_filter_expr() {
            local metric="\$1"
            local operator="\$2"
            local scope="\$3"
            shift 3

            local spec type threshold tag clause pop
            local -a clauses=()
            local -a pop_clauses=()

            for spec in "\$@"; do
                  type="\${spec%%:*}"
                  threshold="\${spec#*:}"

                  enabled "\$threshold" || continue

                  if [[ "\$scope" == "GLOBAL" ]]; then
                        if [[ "\$metric" == "QUAL" ]]; then
                        tag="QUAL"
                        else
                        tag="INFO/\$metric"
                        fi

                        printf -v clause \
                        '(INFO/TYPE="%s" && %s%s%s)' \
                        "\$type" \
                        "\$tag" \
                        "\$operator" \
                        "\$threshold"
                  else
                        pop_clauses=()

                        for pop in "\${POPS[@]}"; do
                        pop_clauses+=(
                              "INFO/\${metric}_\${pop}\${operator}\${threshold}"
                        )
                        done

                        clause=\$(join_by "\$population_join" "\${pop_clauses[@]}")

                        if [[ "${params.vcf_population_fail_mode}" == "ANY" ]]; then
                        clause="(\$clause)"
                        fi

                        printf -v clause \
                        '(INFO/TYPE="%s" && %s)' \
                        "\$type" \
                        "\$clause"
                  fi

                  clauses+=("\$clause")
            done

            if (( \${#clauses[@]} == 0 )); then
                  printf '0'
            else
                  join_by ' || ' "\${clauses[@]}"
            fi
      }

      QUAL_EXPR=\$(make_filter_expr QUAL "<" GLOBAL \
            "SNP:${params.vcf_qual_global_snp}" \
            "INDEL:${params.vcf_qual_global_indel}" \
            "REF:${params.vcf_qual_global_invariant}")

      DP_MIN_EXPR=\$(make_filter_expr DP "<" GLOBAL \
            "SNP:${params.vcf_dp_min_global_snp}" \
            "INDEL:${params.vcf_dp_min_global_indel}" \
            "REF:${params.vcf_dp_min_global_invariant}")

      EH_EXPR=\$(make_filter_expr ExcHet "<" GLOBAL \
            "SNP:${params.vcf_eh_global_snp}" \
            "INDEL:${params.vcf_eh_global_indel}" \
            "REF:${params.vcf_eh_global_invariant}")

      HWE_EXPR=\$(make_filter_expr HWE "<" GLOBAL \
            "SNP:${params.vcf_hwe_global_snp}" \
            "INDEL:${params.vcf_hwe_global_indel}" \
            "REF:${params.vcf_hwe_global_invariant}")

      MAF_EXPR=\$(make_filter_expr MAF "<" GLOBAL \
            "SNP:${params.vcf_maf_global_snp}" \
            "INDEL:${params.vcf_maf_global_indel}" \
            "REF:${params.vcf_maf_global_invariant}")

      NS_EXPR=\$(make_filter_expr NS "<" GLOBAL \
            "SNP:${params.vcf_min_samples_global_snp}" \
            "INDEL:${params.vcf_min_samples_global_indel}" \
            "REF:${params.vcf_min_samples_global_invariant}")

      CR_EXPR=\$(make_filter_expr CR "<" GLOBAL \
            "SNP:${params.vcf_min_callrate_global_snp}" \
            "INDEL:${params.vcf_min_callrate_global_indel}" \
            "REF:${params.vcf_min_callrate_global_invariant}")

      POP_EH_EXPR=\$(make_filter_expr ExcHet "<" POP \
            "SNP:${params.vcf_eh_pop_snp}" \
            "INDEL:${params.vcf_eh_pop_indel}" \
            "REF:${params.vcf_eh_pop_invariant}")

      POP_HWE_EXPR=\$(make_filter_expr HWE "<" POP \
            "SNP:${params.vcf_hwe_pop_snp}" \
            "INDEL:${params.vcf_hwe_pop_indel}" \
            "REF:${params.vcf_hwe_pop_invariant}")

      POP_MAF_EXPR=\$(make_filter_expr MAF "<" POP \
            "SNP:${params.vcf_maf_pop_snp}" \
            "INDEL:${params.vcf_maf_pop_indel}" \
            "REF:${params.vcf_maf_pop_invariant}")

      POP_NS_EXPR=\$(make_filter_expr NS "<" POP \
            "SNP:${params.vcf_min_samples_pop_snp}" \
            "INDEL:${params.vcf_min_samples_pop_indel}" \
            "REF:${params.vcf_min_samples_pop_invariant}")

      POP_CR_EXPR=\$(make_filter_expr CR "<" POP \
            "SNP:${params.vcf_min_callrate_pop_snp}" \
            "INDEL:${params.vcf_min_callrate_pop_indel}" \
            "REF:${params.vcf_min_callrate_pop_invariant}")

      ###########################################
      # Apply soft-filters to VCF (this doesnt yet remove positions)
      ###########################################

      set +e
      bcftools view --threads ${task.cpus} -S ${interval_id}.samples.txt -m2 -M2 -Ou "${variants}" \\
      | bcftools +setGT -Ou -- \\
      -t q \\
      -n . \\
      -i "FORMAT/GQ<${params.vcf_genotype_qual} | FORMAT/DP<${params.vcf_genotype_dp_min} | FORMAT/DP>${params.vcf_genotype_dp_max}" \\
      | bcftools +fill-tags -Ou - -- \\
      -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING' \
      | bcftools +fill-tags -Ou - -- \\
      -S sample_groups.tsv \\
      -t 'NS,MAF,HWE,ExcHet,CR:1=1-F_MISSING' \\
      | bcftools filter -Ou --SnpGap "${params.vcf_dist_indel_global_snp}"  --IndelGap "${params.vcf_dist_indel_global_indel}" \\
      | apply_mask_filters 0  \\
      | bcftools filter -Ou -s QUAL_FAIL -m+ -e "\$QUAL_EXPR" \\
      | bcftools filter -Ou -s DP_MIN_FAIL -m+ -e "\$DP_MIN_EXPR" \\
      | bcftools filter -Ou -s EH_FAIL -m+ -e "\$EH_EXPR" \\
      | bcftools filter -Ou -s HWE_FAIL -m+ -e "\$HWE_EXPR" \\
      | bcftools filter -Ou -s MAF_FAIL -m+ -e "\$MAF_EXPR" \\
      | bcftools filter -Ou -s NS_FAIL -m+ -e "\$NS_EXPR" \\
      | bcftools filter -Ou -s CR_FAIL -m+ -e "\$CR_EXPR" \\
      | bcftools filter -Ou -s POP_EH_FAIL  -m+ -e "\$POP_EH_EXPR" \\
      | bcftools filter -Ou -s POP_HWE_FAIL -m+ -e "\$POP_HWE_EXPR" \\
      | bcftools filter -Ou -s POP_MAF_FAIL -m+ -e "\$POP_MAF_EXPR" \\
      | bcftools filter -Ou -s POP_NS_FAIL  -m+ -e "\$POP_NS_EXPR" \\
      | bcftools filter -Ou -s POP_CR_FAIL  -m+ -e "\$POP_CR_EXPR" \\
      | bcftools view --threads ${task.cpus} -Ob -o tmp.bcf

      # Catch error codes from piped tools so nextflow can retry
      st=("\${PIPESTATUS[@]}")
      set -e
      check_pipeline "\${st[@]}" || exit \$?

      ###########################################
      # Write records passing all soft filters.
      ###########################################

      bcftools view \
            --threads ${task.cpus} \
            --apply-filters PASS \
            --output-type u \
            "tmp.bcf" |
      bcftools annotate \
            --remove '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' \
            --output-type b \
            --output "${interval_id}.filt.bcf"

      bcftools index --threads ${task.cpus} "${interval_id}.filt.bcf"

      # Split the final PASS-only, annotation-cleaned BCF.
      filtered_bcf="${interval_id}.filt.bcf"

      bcftools view \
            --threads ${task.cpus} \
            -v snps \
            -Ob \
            -o "${interval_id}.snp.bcf" \
            "\${filtered_bcf}"

      bcftools view \
            --threads ${task.cpus} \
            -v indels \
            -Ob \
            -o "${interval_id}.indel.bcf" \
            "\${filtered_bcf}"

      bcftools view \
            --threads ${task.cpus} \
            -i 'TYPE="ref"' \
            -Ob \
            -o "${interval_id}.invariant.bcf" \
            "\${filtered_bcf}"

      for type in snp indel invariant; do
            bcftools index \
                  --threads ${task.cpus} \
                  "${interval_id}.\${type}.bcf"
      done

      ###########################################
      # Create QC summaries
      ###########################################

      # Create filter histograms
      create_vcf_filter_histograms \
            "tmp.bcf" \
            "${interval_id}.filter_hist.tsv" \
            "${interval_id}.filter_summary.tsv" \
            100 10 10 \
            "${interval_id}.mask_hist.tsv"
      """
}