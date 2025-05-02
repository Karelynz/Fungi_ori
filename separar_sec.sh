#!/bin/bash

# 🔹 RUTAS DE TRABAJO
INPUT_DIR="./secuencias"
OUTPUT_DIR="./procesadas"
QUALITY_DIR="./fastqc"
TABLE_FILE="conteo_secuencias.csv"

# 🔹 CREAR DIRECTORIOS DE SALIDA
mkdir -p "$OUTPUT_DIR" "$QUALITY_DIR"

# 🔹 ADAPTADORES ITS
ITS5="GGAAGTAAAAGTCGTAACAAGG"
ITS2="GCATCGATGAAGAACGCAGC"
ITS4="TCCTCCGCTTATTGATATGC"

ITS5_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$ITS5"
ITS2_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$ITS2"
ITS4_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$ITS4"

# 🔹 ADAPTADORES CIANOBACTERIAS
CYA106Fw="CGGACGGGTGAGTAACGCGTGA"
CYA781Ra="GACTACTGGGGTATCTAATCCCATT"
CYA781Rb="GACTACAGGGGTATCTAATCCCTTT"

CYA106Fw_FULL="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG$CYA106Fw"
CYA781Ra_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Ra"
CYA781Rb_FULL="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG$CYA781Rb"

# 🔹 ENCABEZADO DE LA TABLA
echo "Muestra,Originales,Hongos,Cianobacterias" > "$TABLE_FILE"

# 🔹 PROCESAR CADA MUESTRA
for file in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    sample_name=$(basename "$file" _R1_001.fastq.gz)
    echo "🔹 Procesando $sample_name..."

    # Verificar existencia de R2
    if [[ ! -f "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" ]]; then
        echo "⚠️  Falta el archivo R2 para $sample_name. Saltando."
        continue
    fi

    # 🔹 Conteo de secuencias originales
    original_count=$(zcat "$file" | echo $((`wc -l`/4)))

    # 🔹 ITS (Hongos)
    cutadapt \
        -g ^$ITS5 -g ^$ITS2 -g ^$ITS4 \
        -g ^$ITS5_FULL -g ^$ITS2_FULL -g ^$ITS4_FULL \
        -o "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" \
        -p "$OUTPUT_DIR/${sample_name}_hongos_R2.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R1_001.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" \
        --discard-untrimmed
    hongos_count=$(zcat "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" | echo $((`wc -l`/4)))

    # 🔹 Cianobacterias
    cutadapt \
        -g ^$CYA106Fw -g ^$CYA781Ra -g ^$CYA781Rb \
        -g ^$CYA106Fw_FULL -g ^$CYA781Ra_FULL -g ^$CYA781Rb_FULL \
        -o "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz" \
        -p "$OUTPUT_DIR/${sample_name}_ciano_R2.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R1_001.fastq.gz" \
        "$INPUT_DIR/${sample_name}_R2_001.fastq.gz" \
        --discard-untrimmed
    ciano_count=$(zcat "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz" | echo $((`wc -l`/4)))

    # 🔹 FASTQC
    fastqc -o "$QUALITY_DIR" "$OUTPUT_DIR/${sample_name}_hongos_R1.fastq.gz" "$OUTPUT_DIR/${sample_name}_ciano_R1.fastq.gz"

    # 🔹 Registrar en la tabla
    echo "$sample_name,$original_count,$hongos_count,$ciano_count" >> "$TABLE_FILE"
done

# 🔹 Generar gráfico en Python
python3 <<EOF
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("$TABLE_FILE")
plt.figure(figsize=(10,6))
plt.bar(df['Muestra'], df['Originales'], label='Originales', alpha=0.6)
plt.bar(df['Muestra'], df['Hongos'], label='Hongos', alpha=0.6)
plt.bar(df['Muestra'], df['Cianobacterias'], label='Cianobacterias', alpha=0.6)
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.ylabel('Número de secuencias')
plt.title('Conteo de secuencias por tipo')
plt.tight_layout()
plt.savefig("conteo_secuencias.png")
print("✅ Gráfico guardado como 'conteo_secuencias.png'")
EOF

echo "✅ Proceso completado. Verifica $OUTPUT_DIR, $QUALITY_DIR y conteo_secuencias.png"