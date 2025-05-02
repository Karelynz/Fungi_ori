#!/bin/bash

# 🔹 Rutas
INPUT_DIR="./secuencias"
OUTPUT_DIR="./procesadas"
mkdir -p "$OUTPUT_DIR" hongos cianos

# 🔹 Adaptadores ITS con overhangs Illumina
ITS5="GGAAGTAAAAGTCGTAACAAGG"
ITS5_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$ITS5"

ITS4="TCCTCCGCTTATTGATATGC"
ITS4_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$ITS4"

# 🔹 Adaptadores para Cianobacterias con overhangs Illumina
CYA106Fw="CGGACGGGTGAGTAACGCGTGA"
CYA781Ra="GACTACTGGGGTATCTAATCCCATT"
CYA781Rb="GACTACAGGGGTATCTAATCCCTTT"

CYA106Fw_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$CYA106Fw"
CYA781Ra_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Ra"
CYA781Rb_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Rb"

# 🔹 Procesamiento de muestras
for file in "$INPUT_DIR"/*_L001_R1_001.fastq.gz; do
  sample_name=$(basename "$file" _L001_R1_001.fastq.gz)
  echo "🔹 Procesando muestra: $sample_name"

  # ITS (Hongos)
  cutadapt \
    -g ^"$ITS5_FULL" \
    -o "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" \
    -p "$OUTPUT_DIR/${sample_name}_hongos_R2.fastq.gz" \
    "$INPUT_DIR/${sample_name}_L001_R1_001.fastq.gz" \
    "$INPUT_DIR/${sample_name}_L001_R2_001.fastq.gz" \
    --discard-untrimmed

  # Cianobacterias
  cutadapt \
    -g ^"$CYA106Fw_FULL" -g ^"$CYA781Ra_FULL" -g ^"$CYA781Rb_FULL" \
    -o "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz" \
    -p "$OUTPUT_DIR/${sample_name}_ciano_R2.fastq.gz" \
    "$INPUT_DIR/${sample_name}_L001_R1_001.fastq.gz" \
    "$INPUT_DIR/${sample_name}_L001_R2_001.fastq.gz" \
    --discard-untrimmed

  # Opcional: análisis de calidad con fastqc
  if command -v fastqc &> /dev/null; then
    fastqc -o "$OUTPUT_DIR" "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz"
    fastqc -o "$OUTPUT_DIR" "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz"
  fi

  # Mover archivos a carpetas específicas
  mv "$OUTPUT_DIR/${sample_name}_hongos_"*.fastq.gz hongos/
  mv "$OUTPUT_DIR/${sample_name}_ciano_"*.fastq.gz cianos/
done

# 🔹 Subir resultados al repositorio
git add hongos/ cianos/
git commit -m "🔼 Secuencias separadas por hongos (ITS5) y cianobacterias (CYA106Fw/CYA781R)"
git push origin main

echo "✅ Proceso completado y archivos subidos a GitHub."