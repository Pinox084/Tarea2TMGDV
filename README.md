Compilar: g++ mlr.cpp -o test
Parametros de ejecución: ./test <n> <eps> <archivo>
Ejemplo de ejecución: ./test 5000000 0.05 Uniform.txt

Compilar codigo de analisis (code.cpp): g++ code.cpp -o code
Parametros de ejecución: ./code <archivo> <eps>
Ejemplo de ejecución: ./code Uniform.txt 0.05 

El archivo code esta dedicado para el análisis extensivo con 10.000 consultas para la función rank, junto con ello obtener los errores y valores exactos de los cuartiles y el rank.

Para ejecución de todos los archivos de texto analizados, usar analisis.sh teniendo ya compilado code.cpp.
Usar chmod +x analisis.sh  y ejecutar.
