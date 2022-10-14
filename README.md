# MTB-HIV
## Overall description

This package includes three partitions, namely the DEA, DCA, and proximity sections. The DEA section is used to generate non-differential interactions (NDIs), while the DCA section is used to generate differential interactions (DIs). The proximity section is applied to drug repurposing based on the gene pair-based distance measures.

## Dependencies

	Perl ( >= 5.0 )
	R  (>= 4.0)
  	Python ( >= 3.8 )

## How to run the DEA and DCA scripts?

Input file

The format of the input file is as follows. You can also check the files in the ./examples/input directory. 

GSE29429_expression.txt  
The first column includes the Entrz IDs of genes and the first row includes the sample IDs.  
ID_REF\tGSM728539\tGSM728540	GSM728541	GSM728542	GSM728543  
1915	  13.19218646	13.27446406	12.9607839	13.01710244	12.95647599  
2597	  11.14440333	11.5689425	10.99561171	11.42964211	11.28005156  
643334	  0	2.125320474	0	0	1.725043085  
9906	  2.729792022	0.486868801	3.523411406	0	0.892429886  
56940	  8.119445668	8.057396116	8.65603151	8.340178134	8.416501993  

PPI.txt  
The first and second columns are Entrz IDs.  
5701	5710  
9343	23016  
3064	26036  
6927	8202  
1525	144100  

GSE29429_sampleInfor.txt  
The first column includes the sample IDs, and the second column means the disease status.  
GSM728545	HC  
GSM728547	HC  
GSM728548	HC  
GSM728539	HMI  
GSM728540	HMI  
GSM728541	HMI  

## Run DEA

Go to ./DEA/bin directory, you will see a run.pl script, and run the following command.

    perl run.pl -expr input_expression_file -sample input_sampleInformation_file -ppi input_ppi_file -disease input_disease_state -out /the/directory/of/output/
    -expr 	the path of expression matrix, such as ../examples/expression.txt
    -sample	the path of sampleInformation
    -ppi 	the path of protein-protein interations
    -disease	the disease state(HMI/MMI/MHCI)
    -out	the path of results
    
    For example:
    perl run.pl -expr ../examples/input/GSE29429_expression.txt -sample ../examples/input/GSE29429_sampleInfor.txt -ppi ../examples/input/PPI.txt -disease HMI -out ../examples/output/

    if you want to detect the differential expressed genes, you can try the following command. (note that the R package Limma is needed)
    Rscript deg.R ../examples/expression.txt ../examples/sampleInfor.txt ../examples/outfile/

Results

You will obtain two output files in the directory /the/directory/of/output/, such as deg.txt, NDI.txt.  
The first file provides the logFC, P.Value, etc. in different conditions.  
The second file provides non-differential interactions.  

Output examples for the first file  
ID_REF	  logFC	AveExpr	t	P.Value	adj.P.Val	B  
9381	  6.25710165951176	4.71547828087234	11.6816199397735	2.40606495322115e-15	7.3384981073245e-11	24.3819451512404  
164668	  3.64056368982549	4.2102917502766	9.62018351259993	1.41090770704937e-12	2.15163425325028e-08	18.3785723457298  
8318	  5.09288604257452	4.19366280910638	9.24617165139942	4.74534620510925e-12	2.62127965596824e-08	17.226138945067  
8208	  3.76792510422157	3.29614476557447	9.24210040425486	4.80884427126567e-12	2.62127965596824e-08	17.2134948500278  

## Run DCA
Go to ./DCA/bin directory, you will see a DCA.R script and run the following command.

	Rscript DCA.R input_expression_file input_sampleInformation_file input_ppi_file HMI /the/directory/of/output/
    input_expression_file     the path of expression matrix, such as ../examples/expression.txt
    input_sampleInformation_file     the path of sampleInformation
    input_ppi_file        the path of protein-protein interations
    /the/directory/of/output/        the path of results
    
    For example:
    Rscript DCA.R ../examples/input/HMI ../examples/input/HC ../data/PPI.txt HMI ../examples/output/

Results

You will obtain four output files in the directory /the/directory/of/output/, such as pairs.txt, GSE29429_disease_pcc.txt, GSE29429_HC_pcc.txt and DI.txt.   
The first file names every gene pair.  
The second and third files provide the PCC value in different conditions.  
The fourth file provides differential interactions.  

Output examples for the second file  
Id	V1  
1	0.682577111307779  
10	-0.106767135104595  
100	-0.247706250913738  
1000	-0.0793841206750866  
10000	0.162843540004409  

Other resources

The normalized expression matrix used in this study are stored in the data directory.  


## How to run the proximity script?
Input file

The format of the input file is as follows. You can also check the files in the ./examples/input directory. 

degree-bin100.txt  
Each line includes the proteins with similar degree measures.  
1	4565	100463486	285382	221914	26290  
76	114112	6076	399473	83259	399726	  
387509	28423	54209	27170	85479	442523  

PPI.pkl  
This file is a dictionary that contains the shortest path length of any two nodes in a PPI network.

drug_data.txt  
The first column is the drug ID, and the following columns are their corresponding targets.

DB00001		2147  	
DB00002		712	713	714	1956	2209	2212	2214	2215  	
DB00004		3559	3560	3561  

HMI_pair.txt 

The first and second columns are the Entrz IDs of HIV-associated gene pairs.

3185	55850  
7322	9921  
1072	5052  
2237	79968  

## Run proximity
Go to ./proximity/bin directory, you will see a proximity.py script and run the following command.

	python proximity.py parm1 parm2 parm3 parm4 parm5
    parm1: The shortest distance file for any two nodes in a PPI network.
    parm2: The association file for drugs and their corresponding targets.
    parm3: The gene pairs related to disease states.
    parm4: The file composed of proteins with similar degree measures.
    parm5: The location of the output file.
    
    For example:
    Python 	proximity.py ../examples/input/PPI.pkl ../examples/input/drug_data.txt ../examples/input/HMI_pair.txt ../examples/input/degree-bin100.txt ../examples/output/HMI

Results

You will obtain four output files in the directory /the/directory/of/output/, such as permutation.txt, proximity.txt, z_score.txt and a directory named permutation that includes 1000 random repetitions.  
The first file provides the average value and standard deviation of each drug.  
The second file provides the three proximity values of each drug.  
The third file provides three Z-scores corresponding to the three proximity values.  

Output examples for the third file

DB00001	1.3972046900472084	0.614549980950421	2.3265382413327567  
DB00002	-1.7558535140358775	-2.150396208830809	0.15986300844872472  
DB00004	1.481542227297944	0.6119972337012101	2.397356726968668  
DB00005	-0.4834025988688505	-1.9133976609548382	2.4905800416108033  

Other resources 

The disease associated gene pairs used in this study are stored in the data directory.

Help and Support

If you have any questions or suggestions, please contact us by email:  jiangyao@webmail.hzau.edu.cn or zhangjiaxuan@webmail.hzau.edu.cn or liurong116@mail.hzau.edu.cn
This software is free for academic use. For commercial use, please contact with the authors.
