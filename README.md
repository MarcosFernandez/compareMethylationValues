===========================
compareMethylationValues.py
===========================

Compares two bedgraph files. Calculates the number of bed positions which are unique to each file, the positions which are 
shared between two files and its methylation values which can be equal or different.

------------------------------------------
Comparison criteria for methylation values
------------------------------------------

A cytosine can be methylated or unmethylated, unfortunately methylation analysis pipelines does not ouput a binary value. Instead,
we got a bedgraph file with an estimated methylation value between 0 and 1.

To consider if two methylation values are equal the next criteria is applied:

1) If the absolute difference between two methylation values is smaller than 0.1 then both values are considered as equal.

2) If both methylation values are in the same methylation cluster then both are considered as equal.

Methylation clusters are defined as, Unmethylated, Undefined and Methylated. Those methylation values between 0 and 0.3 are considered 
as unmethylated, between 0.3 and 0.7 as undefined and between 0.7 and 1 as methylated.

Methylation Clusters and signal ranges:

                |--UnMethylated--|--Undefined--|--Methylated--|
                0               0.3           0.7             1


---------
Licensing
---------

compareMethylationValues script is licensed under GPL. See LICENSE for more information.

-------------
Documentation
-------------

To run the python script two bedgraph files should be provided. These two files should have four fields:
chromosome name, start position, end position and a methylation value between 0 and 1.

These bedgraph files must be sorted, otherwise the results outputted will be wrong.

A label per each file should be provided, it will be used to build stats and barplots.

Example:

    compareMethylationValues.py --first sample1A.sorted.bed --second sample2B.sorted.bed --first-label sample1A --second-label sample2B --png sample1AvsSample2B.png --csv sample2BvsSample1A.csv

Two output files will be generated. Firstone, a csv file with stats regarding unique positions per file and shared positions with same and different methylation values.
And finally a barplot will be generated as png image file.

You can use -r (relaxed comparison) to relax the window location comparison. Then a reported CG can be compared with a C from other source.
         
    |CG|    |C|
     |       |
    |C|     |CG|

------
Author
------

Marcos Fernandez-Callejo at (CNAG/CRG) Centre Nacional d’Anàlisi Genòmica / Centre de Regulació Genòmica.
marcos.fernandez@cnag.crg.eu

