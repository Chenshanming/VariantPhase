process mkdir_bam {
    input:
    path output_path

    output:
    path "$output_path/$params.common.BAM_DIR_NAME"

    script:
    """
    mkdir -p $output_path/$params.common.BAM_DIR_NAME
    """
}

process mkdir_bam_split_by_chr {
    input:
    path output_path

    output:
    path "$output_path/$params.common.BAM_SPLIT_BY_CHR_DIR_NAME"

    script:
    """
    mkdir -p $output_path/$params.common.BAM_SPLIT_BY_CHR_DIR_NAME
    """
}

process mkdir_vcf_split_by_chr {
    input:
    path output_path

    output:
    path "$output_path/$params.common.VCF_SPLIT_BY_CHR_DIR_NAME"

    script:
    """
    mkdir -p $output_path/$params.common.VCF_SPLIT_BY_CHR_DIR_NAME
    """
}

process mkdir_phased_vcf {
  input:
  path output_path

  output:
  path "$output_path/$params.common.PHASED_VCF_DIR_NAME"

  script:
  """
  mkdir -p $output_path/$params.common.PHASED_VCF_DIR_NAME

  """
}

process mkdir_haplotag_bam {
  input:
  path output_path

  output:
  path "$output_path/$params.common.HAPLOTAG_BAM_DIR_NAME"

  script:
  """
  mkdir -p $output_path/$params.common.HAPLOTAG_BAM_DIR_NAME

  """
}

process align {
  input:
  path 'fastq'
  path bam_path

  output:
  path "$bam_path/$params.common.SAM_FILE_NAME"

  script:
  if ("$params.input.ALIGNMENT_METHOD" == "minimap2")
  """
    minimap2 -ax map-$params.input.DATA_TYPE $params.input.REF_FN fastq* > $bam_path/$params.common.SAM_FILE_NAME
  """
}

process sam_to_bam {
  input:
  path bam_path
  path sam_fn

  output:
  path "$bam_path"

  script:
  """
  samtools view -bS --threads $params.input.THREADS $sam_fn > $bam_path/$params.common.BAM_FILE_NAME
  samtools sort $bam_path/$params.common.BAM_FILE_NAME -o $bam_path/$params.common.SORTED_BAM_FILE_NAME

  samtools index $bam_path/$params.common.SORTED_BAM_FILE_NAME
  """
}

process snv_call {
  if("$params.input.SNV_CALL_METHOD" == "clair3")
    conda "$params.input.CONDA_PREFIX/envs/clair3"
  input:
    val sorted_bam_path
    path output_path
  output:
    path "$output_path/snv_call_results"
  script:
  if("$params.input.SNV_CALL_METHOD" == "clair3")
    """
      mkdir -p $output_path/snv_call_results
      run_clair3.sh \
      --bam_fn=$sorted_bam_path/$params.common.SORTED_BAM_FILE_NAME \
      --ref_fn=$params.input.REF_FN\
      --threads=$params.input.THREADS \
      --platform=$params.input.DATA_TYPE \
      --model_path=$params.input.CONDA_PREFIX/envs/clair3/bin/models/$params.input.DATA_TYPE \
      --output=$output_path/snv_call_results
    """
}

process filterVariant {
  input:
  
  output:
  
  script:
  """
  
  """
} 

process split_by_chr {
  input:
  path vcf_split_by_chr_path
  path snv_results_path
  path bam_split_by_chr_path
  path sorted_bam_path
  each chr_idx
  output:
  val "split.done"
  script:
  if("$params.input.SNV_CALL_METHOD" == "clair3")
  """
  tabix -f -p vcf $snv_results_path/$params.common.CLAIR3_VCF_NAME
  bcftools view -r chr${chr_idx} $snv_results_path/$params.common.CLAIR3_VCF_NAME > $vcf_split_by_chr_path/chr${chr_idx}.vcf

  samtools view -b -h $sorted_bam_path/$params.common.SORTED_BAM_FILE_NAME chr${chr_idx} > $bam_split_by_chr_path/chr${chr_idx}.bam
  samtools index $bam_split_by_chr_path/chr${chr_idx}.bam

  """
}
  
process phase {
  input:
  path phased_vcf_path
  path ref_fn
  path bam_split_by_chr_path
  path vcf_split_by_chr_path
  each chr_idx
  val split_done
  output:
  val "phased.done"
  script:
  """
  samtools faidx $ref_fn
  whatshap phase \
  --output $phased_vcf_path/chr${chr_idx}.vcf \
  --reference $ref_fn \
  --chromosome chr${chr_idx} \
  --distrust-genotypes \
  --ignore-read-groups \
  $vcf_split_by_chr_path/chr${chr_idx}.vcf \
  $bam_split_by_chr_path/chr${chr_idx}.bam

  bgzip -c $phased_vcf_path/chr${chr_idx}.vcf > $phased_vcf_path/chr${chr_idx}.vcf.gz
  tabix -p vcf $phased_vcf_path/chr${chr_idx}.vcf.gz
  """
}

process haplotag {

  input:
  path haplotag_bam_path
  path ref_fn
  path bam_split_by_chr_path
  path phased_vcf_path
  each chr_idx
  val phased_done
  output:

  script:
  """
  samtools faidx $ref_fn
  whatshap haplotag \
  --output $haplotag_bam_path/chr${chr_idx}.bam \
  --reference $ref_fn \
  --ignore-read-groups \
  --regions chr${chr_idx} \
  --output-haplotag-list $haplotag_bam_path/chr${chr_idx}.tsv \
  $phased_vcf_path/chr${chr_idx}.vcf.gz \
  $bam_split_by_chr_path/chr${chr_idx}.bam
  """

}
  
process test {
    conda "$CONDA_PREFIX/envs/clair3"

    input:
    val snv_call_method
    path fastaq
    path test_log
    output:
    
    script:
    """
    echo $fastaq > test.log
    """
}   
    
workflow  {
  FASTQ_PATH=Channel.fromPath("$params.input.FASTQ_PATH/*.fastq")
  REF_FN=Channel.fromPath("${params.input.REF_FN}")
  OUTPUT_PATH=Channel.fromPath("${params.input.OUTPUT_PATH}")

  chr_idx = Channel.fromList([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'])

  // align

  bam_path = mkdir_bam(OUTPUT_PATH)

  //sam_fn = align(FASTQ_PATH, bam_path)

  //sorted_bam_path = sam_to_bam(bam_path, sam_fn)

  // call snv

  snv_results_path = snv_call(bam_path, OUTPUT_PATH)

  // preprocess before phasing

  bam_split_by_chr_path = mkdir_bam_split_by_chr(OUTPUT_PATH)

  vcf_split_by_chr_path = mkdir_vcf_split_by_chr(OUTPUT_PATH)

  split_done = split_by_chr(vcf_split_by_chr_path, snv_results_path, bam_split_by_chr_path, bam_path, chr_idx)

  // phasing
  phased_vcf_path = mkdir_phased_vcf(OUTPUT_PATH)

  phased_done = phase(phased_vcf_path, REF_FN, bam_split_by_chr_path, vcf_split_by_chr_path, chr_idx, split_done)

  // haplotaging
  haplotag_bam_path = mkdir_haplotag_bam(OUTPUT_PATH)

  haplotag(haplotag_bam_path, REF_FN, bam_split_by_chr_path, phased_vcf_path, chr_idx, phased_done)


}   