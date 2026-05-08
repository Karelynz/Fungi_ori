# Fungi_ori: Procesamiento de secuencias fúngicas y cianobacterianas

Este repositorio contiene un pipeline automatizado para el recorte, análisis y visualización de secuencias de ADN dirigidas a hongos (ITS) y cianobacterias (16S), incluyendo control de calidad y generación de reportes.

---

## Estructura del repositorio

- `separar_sec.sh`: script principal que ejecuta todo el flujo.
- `hongos/`: archivos FASTQ recortados para hongos.
- `cianos/`: archivos FASTQ recortados para cianobacterias.
- `logs/`: registros de recorte por muestra y resumen global.
- `multiqc/`: reporte interactivo de calidad.
- `resumen_cutadapt.png`: gráfico resumen de recorte de adaptadores.
- `multiqc_report.html`: reporte MultiQC visual.

---

##  Ejecución del pipeline

```bash
chmod +x separar_sec.sh
./separar_sec.sh
Esto ejecuta:

Recorte de adaptadores con cutadapt
Análisis de calidad con fastqc (opcional) y multiqc (si está disponible)
Resumen y gráficos de recorte por muestra

 Resumen del recorte de adaptadores

 Informe de calidad (FastQC + MultiQC)

Ver reporte completo

 Requisitos

cutadapt
fastqc (opcional pero recomendado)
multiqc (reporte agregado de calidad)
python3, matplotlib, pandas (para generar el gráfico)
Instalación rápida:

pip install cutadapt fastqc multiqc matplotlib pandas
 Adaptadores utilizados

Hongos:
ITS5 (R1): GGAAGTAAAAGTCGTAACAAGG
ITS4 (R2): TCCTCCGCTTATTGATATGC
Cianobacterias:
CYA106Fw (R1): CGGACGGGTGAGTAACGCGTGA
CYA781Ra (R2): GACTACTGGGGTATCTAATCCCATT
CYA781Rb (R2): GACTACAGGGGTATCTAATCCCTTT
Autora

Karen Nunez (github.com/Karelynz)
Reserva de la Biosfera de Mapimí – Análisis de comunidades microbianas en biocostras
Desarrollado durante mi tesis doctoral - IPICyT
