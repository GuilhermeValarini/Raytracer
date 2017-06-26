Neste trabalho temos 3 implementações de um algoritmo de Raytracing: Serial, paralello por OpenMP e paralelo por CUDA

Cada arquivo de entrada é uam cena a ser renderizada definida no seguinte formato:
linhas: Valores no arquivo                            -> Explicação
1     : w h                                           -> width e height da imagem a ser renderizada
2     : s l                                           -> Número de esferas e de luzes presentes na cena
                                                         (Obs: em nossos ambiente temos apenas 1 luz)
s     : x y z raio red green blue reflection          -> Coordenadas (x,y,z) e o Raio da esfera. Cor da esfera (r,g,b) e
                                                         um indice de reflectância da esfera
l     : x y z raio red green blue reflection emission -> Coordenadas (x,y,z) e o Raio da luz. Cor da luz (r,g,b) e um indice
                                                         de reflectância da luz. Por último temos a intensidade da luz

                                                         O nome destes arquivos serão os mesmos dos arquivos de entrada que os geraram.
Cada entrada gera uma imagem de saída que estará presente dentro da pasta de cada implementação.
