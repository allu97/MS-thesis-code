% cd 'C:\Users\Allu\Documents\GitHub\XRDidentifier'
test = readtable('testing_curve.csv');
T = readtable('learning_curve.csv');
figure
hold on
y = test.acc_test > 0;
plot(T.Var1,T.acc_train,'b','linewidth',2)
plot(T.Var1(y),T.acc_val(y),'r','linewidth',2)
plot(test.Var1(y),test.acc_test(y),'g','linewidth',2)
plot(test.Var1(y),T.acc_val(y),'ko','linewidth',1.5)
legend('Train','Validation','Testing','Model saved','location','se')
ylabel('Accuracy (%)')
xlabel('Epoch')
title('CNN learning curve')
grid on

figure;hold on
plot(T.Var1,T.loss_train,'b')
plot(T.Var1,T.loss_val,'r')
