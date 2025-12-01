#!/bin/bash

# ============================================================
#  Script para ejecutar el análisis experimental del sketch MRL
#  para los cuatro archivos entregados, con epsilon = 0.1 y 0.05
# ============================================================

# Archivos a procesar
FILES=("chicago2015.txt" "chicago2016.txt" "Uniform.txt" "Log-normal.txt")

# Valores de epsilon
EPSILONS=("0.1" "0.05")


EXEC="./code"


if [ ! -f "$EXEC" ]; then
    echo "ERROR: No se encontró el ejecutable '$EXEC'"
    echo "Compila con: g++ code.cpp -o code"
    exit 1
fi

mkdir -p resultados

echo "==============================================="
echo "   EJECUTANDO EXPERIMENTOS MRL AUTOMÁTICOS"
echo "==============================================="
echo ""

for file in "${FILES[@]}"; do
    for eps in "${EPSILONS[@]}"; do

        OUT="resultados/${file}_${eps}.txt"

        echo "Ejecutando análisis para archivo='$file' con epsilon=$eps ..."
        
        
        $EXEC "$file" "$eps" > "$OUT"

        echo " → Guardado en: $OUT"
        echo ""

    done
done

echo "==============================================="
echo "  Todos los experimentos han sido completados."
echo "  Revisa la carpeta 'resultados/'."
echo "==============================================="
