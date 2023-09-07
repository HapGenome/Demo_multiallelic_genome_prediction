# multi_allelic_genome_prediction
This project is designed to demo how to perform multi-allelic genome prediction, which includes the following steps:
  Step1: Pull/Get the multi-alleic data, eg. ancestral haplotype from CSW
  Step2: Perform numeric coding for this data. A R function called getDummyMatrixForPbhData.R is a example of how we encode the ancestral haplotypes. And example is shared (oneHotDataPred_20230829_R.ipynb)
  Step3: Run genome prediction models. A example of using Elastic net, LightGBM and DNN has been provided (hap_inf_pred_2020.ipynb)
