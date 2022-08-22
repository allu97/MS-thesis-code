# MS-thesis-code
My master's thesis code for data handling (NaN removal, removing variables, removing data rows/objects, visualization, testing different machine learning methods for classification, different pre-treatment, and so on.

The data is divided into two parts: XRF (X-ray fluorescence) and FTIR (Fourier-transform infrared spectroscopy).

Thesis (paper): https://urn.fi/URN:NBN:fi-fe2022060242271
Introduction video: https://youtu.be/1JhMyr8wmUo

# Citation
The CNN with metric learning used for FTIR data is based on the work fetched from https://github.com/ma921/XRDidentifier and based on the article from https://towardsdatascience.com/automatic-spectral-identification-using-deep-metric-learning-with-1d-regnet-and-adacos-8b7fb36f2d5f 

My work for the 1-D CNN with metric learning is modified to fetch data from .csv files and use top-1 classification instead of top-5. My version does not use data augmentation since there are many samples for each class.

Other work:

Papers

AdaCos: https://arxiv.org/abs/1905.00292
1D-RegNet: https://arxiv.org/abs/2008.04063
Physics-informed data augmentation: https://arxiv.org/abs/1811.08425v2
Sparsely-gated layer: https://arxiv.org/abs/1701.06538

Implementation

AdaCos: https://github.com/4uiiurz1/pytorch-adacos/blob/master/metrics.py
1D-RegNet: https://github.com/hsd1503/resnet1d
Physics-informed data augmentation: https://github.com/PV-Lab/autoXRD
Top k accuracy: https://gist.github.com/weiaicunzai/2a5ae6eac6712c70bde0630f3e76b77b
Angular Penalty Softmax Loss: https://github.com/cvqluu/Angular-Penalty-Softmax-Losses-Pytorch
Sparsely-gated layer: https://github.com/davidmrau/mixture-of-experts
