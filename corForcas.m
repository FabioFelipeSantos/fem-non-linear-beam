function Cor = corForcas(fi)
Cores = [255, 0, 0;
         0, 179, 0;
         255, 0, 255;
         153, 0, 204;
         255, 153, 51;
         204, 51, 0;
         0, 51, 102;
         204, 153, 0;
         204, 0, 102;
         0, 153, 255] / 255;
Cores = [Cores;Cores;Cores;Cores;Cores];
Cor = Cores(fi, :);