<h2>Installation</h2>
<p><b>Library and Django Version is super important !</b></p>
<p><b>Pre-STEP1</b>: install Django(v1.11.4), and MySQL</p>
<p><b>Pre-STEP2</b>: pip install numpy scipy pandas matplotlib matplotlib_venn networkx gseapy lifelines==0.11.2 celery seaborn</p>

<p><b>STEP1</b>: Please use git download for web application, and for database, you have to get MySQL dump file and PCTA datafile from the official sites</p>
Link : http://www.thepcta.org
<p><b>STEP2</b>: Create PCTA database on MySQL, database name should be 'pcta' and import sql file on your MySQL.
(command mysql -uroot -p pcta < pcta.sql)</p>
<p><b>STEP3</b>: Move pcta_expression_data.csv to 'user_data' folder in your downloaded PCTA django folder.</p>
<p><b>STEP4</b>: Change 'your_googlegroup' to your google group name in 'qna.html' file</p>
<p>document.getElementById('forum_embed').src = 'https://groups.google.com/forum/embed/?place=forum/your_googlegroup'<p>
<p><b>STEP5</b>: run python manage.py migrate</p>
<p><b>STEP6</b>: run python manage.py runserver localhost:8000</p>
<p><b>STEP7</b>: run celery -A pcta_updated worker</p>

<p><b>Additional STEP</b>: Install Nginx, uWSGI, and please use socket file and uWSGI init file in uWSGI folder to link uWSGI and Nginx.</p>

<h2>Manual</h2>
        <h4>Type of analysis</h4>
        <p><b>Association analysis</b> : This analysis is composed of 3 major plots. Waterfall plot, Violin plot and Histogram. Each plots are made by PCTA dataset, and you can see expression trends of your gene or gene list.</p>
        <p><b>Correlation analysis</b> : This analysis is composed of scatter plots and regression line. As a result, you can check correlation between 2 sets of your input through PCTA dataset. Correlation statistic is spearman rank sum.</p>
        <p><b>Set analysis</b> : This analysis is a big category of 2 different analyses. Gene Set Enrichment Analysis and Master Regulator Analysis. You will get GSEA result of your gene set and its master regulator candidates through PCTA dataset.</p>

<h4>Input</h4>
        <p>Entrez ID and official gene symbol in this version.</p>
        <p><b>Pathway input</b> : Click Pathway input button above input box, then you can see dialog box to choose pathway. Choose one of them, and Click enter button</p>
        <p><b>Association analysis</b> : Copy and paste a gene and gene set. If input is gene set, it will be calculated to Z score.</p>
        <p><b>Correlation analysis</b> : Copy and paste a gene and gene set in 2 different input boxes. If input is gene set, it will be calculated to Z score. Additionally, Input name can be customized</p>
        <p><b>Set analysis</b> : Copy and paste gene set only and input name can be customized. Set analysis needs more than 5 genes for the input to increase accuracy.</p>
        <p><b>Supplemnet</b> : If you enter a set of genes, it will be calculated Z score(Gene set Z score) automatically not original expression values.</p>

<h4>Option</h4>
        <p><b>Association analysis</b> : This analysis has 4 options. Disease course, PCS, PAM50 and BCR. Disease course option will divide PCTA samples by Gleason score (GS<7, GS=7, GS>7, mCRPC). PCS is one of prostate cancer stratification system, and it will categorized samples by PCS1, PCS2 and PCS3. PAM50 is one of prostate cancer stratification system as similar as PCS, and it will divide PCTA dataset by Luminal A, Luminal B and Basal. BCR means Biochemical Recurrence Free analysis, and you can get Kaplan-Meier plot and Cox Proportional Hazard Analysis for your input.</p>
        <p><b>Correlation analysis</b> : This analysis has 3 options. Disease course, PCS, PAM50 and BCR. Disease course option will divide PCTA samples by Gleason score (GS<7, GS=7, GS>7, mCRPC). PCS is one of prostate cancer stratification system, and it will categorized samples by PCS1, PCS2 and PCS3. PAM50 is one of prostate cancer stratification system as similar as PCS, and it will divide PCTA dataset by Luminal A, Luminal B and Basal.</p>
        <p><b>Set analysis</b> : This analysis has bi-sampling option(ex. GS<7 versus Others in Disease Course). Bi-sampling option is 3 major categories, Disease course, PCS and PAM50.</p>

<h4>Result</h4>
You can get results of every analysis by clicking Image Download  and Table Download at the top of result screen.

<h4>Download raw data</h4>
You can download PCTA dataset and its clinical data in Download section.

<h4>Reference</h4>

Main spec and library information
Main workframe : Django
1. GSEA : gseapy for python (<a href="http://gseapy.readthedocs.io/en/latest/">Link</a>)
2. Network plot : networkx for python (<a href="https://networkx.github.io/">Link</a>)
3. Survival analysis : lifelines for python (<a href="http://lifelines.readthedocs.io/en/latest/#">Link</a>)
4. Message passing library : celery (<a href="http://www.celeryproject.org/">Link</a>))
5. ETC library : pandas, scipy, numpy, seaborn, matplotlib
6. Database : MySQL


