data = [3015.04
104.647
65.5927
16.7684
9.92434
6.3368
4.10323
2.50191
2.1816
1.91978
1.38394
1.15388
1.02549
0.994836
0.871269
0.819393
0.742435
0.685258
0.627377
0.607419
0.585025
0.519298
0.492307
0.459208
0.423575
0.402649
0.36863
0.342275
0.328659
0.314581
0.282152
0.269562
0.265521
0.255774
0.232769
0.224623
0.208705
0.173441
0.172148
0.159333
0.144937
0.139018
0.1335
0.12612
0.111007
0.109431
0.100254
0.0932766
0.0865681
0.0682743
0.0613558
0.0563453
0.0537711
0.0500841
0.0380819
1.94E-07
];

%Pre allocate to save time
cont = zeros(56,1);
%calculate cumulative sums up to i
for i=1:56
    for j=1:i
        cont(i) = cont(i) + data(j);
    end
end

%find total so that the percentage of total can be found
tot = 0;
for i=1:56
    tot = tot + data(i);
end


for i =1:56
    cont(i) = cont(i)/tot;
end

plot(cont);
title("Explained variance ratio");
xlabel("Number of Eigenvalues");
ylabel("Explained Variance");
