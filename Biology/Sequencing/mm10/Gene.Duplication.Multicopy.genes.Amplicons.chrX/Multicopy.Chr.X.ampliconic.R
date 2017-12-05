
setup_MarkdownReports(OutDir = '/Users/abelvertesy/Github_repos/TheCorvinas/Biology/Sequencing/mm10/Gene.Duplication.Multicopy.genes.Amplicons.chrX/')


Multicopy.GeneFamilies.chrX.mm10 = c( 'EG668965', 'Gmcl1', 'Ssxb', 'Fthl17', 'Zfp161', 'Slx', 'E330016L19Rik', 'Rhox', 'Cxx1', '4930527E24Rik', 'Magea', '4930567H17Rik', 'Xlr', 'Sstx', 'EG434797', 'Obp1', 'Mageb', 'Zxd', 'Dmrtc1b', 'Pabpc1l2', 'Tgif2lx', 'Srsx', 'LOC665542', 'Pramel3', 'Ott', 'LOC207318', '4921511M17Rik', 'Magea')
Multicopy.Genes.chrX.mm10.NonAmpliconic = c( '1700012L04Rik', '1700020N15Rik', '1700042B14Rik', '1700003E24Rik', '1700010D01Rik', 'MGC58426', '4933434C23Rik', '4930524N10Rik', 'LOC278181', 'LOC245376')
Multicopy.Genes.chrX.mm10.Ampliconic = c( 'Gmcl1l', 'EG668965', 'Ssxb', 'Fthl17', 'Zfp161', 'Slx', 'E330016L19Rik', 'Rhox', '4930527E24Rik', '4930567H17Rik', 'Magea', 'Zxd', 'Dmrtc1b', 'Pabpc1l2', 'Tgif2lx', 'Srsx', 'LOC665542', 'Pramel3', 'Ott', 'LOC207318', '4921511M17Rik', 'Mageb', 'EG434797', 'Xlr')

l(Multicopy.GeneFamilies.chrX.mm10)
l(Multicopy.Genes.chrX.mm10.NonAmpliconic)
l(Multicopy.Genes.chrX.mm10.Ampliconic)


Multicopy.Genes.chrX.mm10.Merged = c(Multicopy.Genes.chrX.mm10.NonAmpliconic, Multicopy.Genes.chrX.mm10.Ampliconic)

UNIONNN = union(Multicopy.GeneFamilies.chrX.mm10,Multicopy.Genes.chrX.mm10.Merged)
toClipboard(UNIONNN)

# YAY done
Updated.GeneIDs.for.all.multicopy.genes.chrX = c( 'Gm9427', 'Gmcl1', 'Ssxb', 'Fthl17', 'Zfp161', 'Slx', 'E330016L19Rik', 'Rhox13', 'Rhox2a', 'Rhox2b', 'Rhox2c', 'Rhox2d', 'Rhox2e', 'Rhox2f', 'Rhox2g', 'Rtl8a', 'Rtl8b', 'Rtl8c', 'Slxl1', 'Magea10', 'Magea2', 'Magea3', 'Magea5', 'Magea8', 'Magea4', 'Magea6', 'Gm14725', 'Xlr', 'Sstx', 'Gm14744', 'Mageb', 'Zxdc', 'Dmrtc1b', 'Pabpc1l2', 'Tgif2lx1', 'Gm17503')

