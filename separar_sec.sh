#!/bin/bash

# ============================================
# SEPARACION DE SECUENCIAS Y SUBIDA A GITHUB
# ============================================

# 🔹 Rutas base
INPUT_DIR="./secuencias"
OUTPUT_DIR="./procesadas"
QUALITY_DIR="./calidad"
GIT_REPO_DIR="./Fungi_ori"
HONGOS_DIR="$GIT_REPO_DIR/hongos"
CIANO_DIR="$GIT_REPO_DIR/cianos"

# 🔹 Crear carpetas necesarias
mkdir -p "$OUTPUT_DIR" "$QUALITY_DIR" "$HONGOS_DIR" "$CIANO_DIR"

# 🔹 Adaptadores para ITS (Hongos)
ITS1="CTTGGTCATTTAGAGGAAGTAA"
ITS2="GCTGCGTTCTTCATCGATGC"
ITS3="GCATCGATGAAGAACGCAGC"

ITS1_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$ITS1"
ITS2_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$ITS2"
ITS3_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$ITS3"

# 🔹 Adaptadores para Cianobacterias
CYA106Fw="CGGACGGGTGAGTAACGCGTGA"
CYA781Ra="GACTACTGGGGTATCTAATCCCATT"
CYA781Rb="GACTACAGGGGTATCTAATCCCTTT"

CYA106Fw_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$CYA106Fw"
CYA781Ra_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Ra"
CYA781Rb_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Rb"

# 🔹 Recorrer muestras
for file in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    sample_name=$(basename "$file" _R1_001.fastq.gz)
    echo "🔹 Procesando $sample_name..."

    # Verificar existencia de pares
    if [[ ! -f "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" ]]; then
        echo "⚠️ Pares no encontrados para $sample_name."
        continue
    fi

    # 🔸 Separar secuencias de HONGOS
    cutadapt \
        -g ^$ITS1_FULL \
        -g ^$ITS2_FULL \
        -g ^$ITS3_FULL \
        -o "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" \
        -p "$OUTPUT_DIR/${sample_name}_hongos_R2.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R1_001.fastq.gz" "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" \
        --discard-untrimmed

    # 🔸 Separar secuencias de CIANOBACTERIAS
    cutadapt \
        -g ^$CYA106Fw_FULL \
        -g ^$CYA781Ra_FULL \
        -g ^$CYA781Rb_FULL \
        -o "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz" \
        -p "$OUTPUT_DIR/${sample_name}_ciano_R2.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R1_001.fastq.gz" "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" \
        --discard-untrimmed

    # 🔸 Análisis de calidad posterior (opcional si fastqc está disponible)
    if command -v fastqc &> /dev/null; then
        fastqc -o "$QUALITY_DIR" "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz"
    fi

    # 🔸 Mover archivos separados al repositorio Git
    mv "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" "$HONGOS_DIR/"
    mv "$OUTPUT_DIR/${sample_name}_hongos_R2.fastq.gz" "$HONGOS_DIR/"
    mv "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz" "$CIANO_DIR/"
    mv "$OUTPUT_DIR/${sample_name}_ciano_R2.fastq.gz" "$CIANO_DIR/"

done

# 🔹 Subir a GitHub
cd "$GIT_REPO_DIR"
git pull origin main
git add hongos/ cianos/
git commit -m "🔹 Añadidas secuencias separadas por tipo (hongos y cianos)"
git push origin main

echo "✅ Secuencias separadas y subidas a GitHub correctamente."