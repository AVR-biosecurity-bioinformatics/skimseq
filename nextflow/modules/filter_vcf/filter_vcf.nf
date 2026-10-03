process FILTER_VCF {

    tag "${interval_hash}"
    conda "${moduleDir}/environment.yml"

      input:
      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path(bcf),
            path(csi)

      path mask_bed
      path popmap
      path missing_summary

      output:
      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path("${interval_hash}.filt.bcf"),
            path("${interval_hash}.filt.bcf.csi"),
            emit: all_bcf

      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path("${interval_hash}.snp.bcf"),
            path("${interval_hash}.snp.bcf.csi"),
            emit: snp_bcf

      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path("${interval_hash}.indel.bcf"),
            path("${interval_hash}.indel.bcf.csi"),
            emit: indel_bcf

      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path("${interval_hash}.invariant.bcf"),
            path("${interval_hash}.invariant.bcf.csi"),
            emit: invariant_bcf

      tuple val(interval_hash),
            path(interval_bed),
            path(bed_tbi),
            path("${interval_hash}.metrics.tsv.gz"),
            emit: metrics

      path("${interval_hash}.samples.txt"),
            emit: samples_to_keep

      script:
      """
      #!/usr/bin/env bash
      set -euo pipefail
      
      # Source dependent functions
      source "\$(command -v functions.sh)"

      # Make sure mask file is sorted and unique (and 0-based, half-open)
      sort -k1,1 -k2,2n -k3,3n ${mask_bed} | uniq > vcf_masks.bed

      # Find samples above the missing fraction filter
      awk -v thr="${params.sample_max_missing}" 'NR==1 {next} \$4!="NA" && (\$4+0) < thr {print \$1}' "${missing_summary}" > ${interval_hash}.samples.txt

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
      ' "${interval_hash}.samples.txt" "${popmap}" > sample_groups.tsv

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

      set +e
      bcftools view --threads ${task.cpus} -S ${interval_hash}.samples.txt -m2 -M2 -Ou "${bcf}" \
      | bcftools +setGT -Ou -- \
      -t q \
      -n . \
      -i "FORMAT/GQ<${params.vcf_genotype_qual} || FORMAT/DP<${params.vcf_genotype_dp_min} || FORMAT/DP>${params.vcf_genotype_dp_max}" \
      | bcftools +fill-tags -Ou - -- \
      -t 'AC,AN,NS,MAF,F_MISSING,HWE,ExcHet,TYPE,CR:1=1-F_MISSING' \
      | bcftools +fill-tags -Ou - -- \
      -S sample_groups.tsv \
      -t 'NS,MAF,HWE,ExcHet,CR:1=1-F_MISSING' \
      | bcftools filter -Ou --SnpGap "${params.vcf_dist_indel_global_snp}"  --IndelGap "${params.vcf_dist_indel_global_indel}" \
      | bcftools filter -Ou -s MASK_FAIL -m+ -M vcf_masks.bed \
      | bcftools filter -Ou -s QUAL_FAIL -m+ -e "\$QUAL_EXPR" \
      | bcftools filter -Ou -s DP_MIN_FAIL -m+ -e "\$DP_MIN_EXPR" \
      | bcftools filter -Ou -s EH_FAIL -m+ -e "\$EH_EXPR" \
      | bcftools filter -Ou -s HWE_FAIL -m+ -e "\$HWE_EXPR" \
      | bcftools filter -Ou -s MAF_FAIL -m+ -e "\$MAF_EXPR" \
      | bcftools filter -Ou -s NS_FAIL -m+ -e "\$NS_EXPR" \
      | bcftools filter -Ou -s CR_FAIL -m+ -e "\$CR_EXPR" \
      | bcftools filter -Ou -s POP_EH_FAIL  -m+ -e "\$POP_EH_EXPR" \
      | bcftools filter -Ou -s POP_HWE_FAIL -m+ -e "\$POP_HWE_EXPR" \
      | bcftools filter -Ou -s POP_MAF_FAIL -m+ -e "\$POP_MAF_EXPR" \
      | bcftools filter -Ou -s POP_NS_FAIL  -m+ -e "\$POP_NS_EXPR" \
      | bcftools filter -Ou -s POP_CR_FAIL  -m+ -e "\$POP_CR_EXPR" \
      | bcftools view --threads ${task.cpus} -Ob -o tmp.bcf

      # Catch error codes from piped tools so nextflow can retry
      st=("\${PIPESTATUS[@]}")
      set -e
      check_pipeline "\${st[@]}" || exit \$?

      # Write records passing all soft filters.
      bcftools view \
            --threads ${task.cpus} \
            --apply-filters PASS \
            --output-type u \
            "tmp.bcf" |
      bcftools annotate \
            --remove '^INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING,INFO/HWE,INFO/ExcHet,INFO/TYPE,INFO/CR' \
            --output-type b \
            --output "${interval_hash}.filt.bcf"

      bcftools index --threads ${task.cpus} "${interval_hash}.filt.bcf"

      # Split the final PASS-only, annotation-cleaned BCF.
      filtered_bcf="${interval_hash}.filt.bcf"

      bcftools view \
            --threads ${task.cpus} \
            -v snps \
            -Ob \
            -o "${interval_hash}.snp.bcf" \
            "\${filtered_bcf}"

      bcftools view \
            --threads ${task.cpus} \
            -v indels \
            -Ob \
            -o "${interval_hash}.indel.bcf" \
            "\${filtered_bcf}"

      bcftools view \
            --threads ${task.cpus} \
            -i 'TYPE="ref"' \
            -Ob \
            -o "${interval_hash}.invariant.bcf" \
            "\${filtered_bcf}"

      for type in snp indel invariant; do
            bcftools index \
                  --threads ${task.cpus} \
                  "${interval_hash}.\${type}.bcf"
      done

      # Population-specific INFO tags vary with the retained populations, so discover
      # only those tags dynamically. Global tags are written explicitly below.
      mapfile -t POP_INFO_TAGS < <(
      bcftools view --header-only "tmp.bcf" |
      awk '
            /^##INFO=<ID=/ {
                  tag = \$0
                  sub(/^##INFO=<ID=/, "", tag)
                  sub(/,.*/, "", tag)

                  if (tag ~ /^(NS|MAF|HWE|ExcHet|CR)_[^_]+\$/)
                  print tag
            }
      '
      )

      header=\$'CHROM\\tPOS\\tFILTER\\tQUAL\\tTYPE\\tDP\\tExcHet\\tHWE\\tMAF\\tNS\\tCR'
      format='%CHROM\\t%POS\\t%FILTER\\t%QUAL\\t%INFO/TYPE\\t%INFO/DP\\t%INFO/ExcHet\\t%INFO/HWE\\t%INFO/MAF\\t%INFO/NS\\t%INFO/CR'

      for tag in "\${POP_INFO_TAGS[@]}"; do
      header+=\$'\\t'"\$tag"
      format+=\$'\\t'"%INFO/\$tag"
      done

      format+='\\n'

      {
      printf '%s\\n' "\$header"
      bcftools query \
            --format "\$format" \
            "tmp.bcf"
      } |
      bgzip --stdout > "${interval_hash}.metrics.tsv.gz"
      """
}