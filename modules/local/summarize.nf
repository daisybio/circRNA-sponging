process SUMMARIZE {
  input:
    file(circRNA_files)
  
  output:
    file("circRNA_counts.tsv")
  
  script:
  """
    ${projectDir}/scripts/circRNA/summarize.py --input ${circRNA_files} --output circRNA_counts.tsv
  """
}