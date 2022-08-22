[coeff,score,latent,tsquared,explained] = pca(Reflectances);
idx1 = Class==1;
idx2 = Class==2;
idx3 = Class==3;
idx4 = Class==4;

figure
biplot('Scores',score(:,1:2));
figure;hold on
plot(score(idx1,1),score(idx1,2),'r*')
plot(score(idx2,1),score(idx2,2),'gx')
plot(score(idx3,1),score(idx3,2),'bo')
plot(score(idx4,1),score(idx4,2),'k*')

figure;hold on
plot(score(idx1,1),score(idx1,4),'r*')
plot(score(idx2,1),score(idx2,4),'gx')
plot(score(idx3,1),score(idx3,4),'bo')
plot(score(idx4,1),score(idx4,4),'k*')
