import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('learning_curve.csv', index_col=0)
dt = pd.read_csv('testing_curve.csv',index_col=0)
plt.plot(df.index, df['acc_train'], label='train')
plt.plot(df.index, df['acc_val'], label='validation')
y = dt['acc_test']
print(y.shape)
for row in range(y.shape[0]):
    if y[row] > 0:
        plt.plot(row, y[row],'bo')
plt.legend()
plt.xlabel('epoch')
plt.ylabel('accuracy (%)')
plt.savefig('learning_curve.png')