#%%
import math
import matplotlib.pyplot as plt
import numpy as np

rates = np.array((0.705, 0.875, 1.043, 1.235, 1.445))
ttms = np.array((0.5, 1, 1.5, 2, 2.5))

swapRates = np.column_stack((ttms, rates))

# %%
class bootstrap():
    
    def discountCurve(N, Freq, x):
        """
        solve for df in: 
            1 = sum(couponAmount(i)/couponFreq * df(i)) + (1 + couponAmount(n)) * df(n)
        """
        n = len(x)
        output = np.zeros(n)
        for i in range(0, n):
            coup = sum([x[i][1] / Freq * output[t] for t in range(0, i)])
            output[i] = (N - coup) / (N + x[i][1] / Freq)
        return output
    
    def spotRates(N, Freq, x):
        """
        solve for r in:
            1 = (1 / df(i)) ^(1 / ttm * Freq)
        """
        n = len(x)
        output = np.zeros(n)
        for i in range(0, n):
            output[i] = (((N / x[i][1]) ** (1 / (x[i][0] * Freq))) - 1) * Freq
        return output
    
    def forwardRates(N, Freq, x):
        """
        solve for FR in:
            (1 + r(i) / Freq) ^ (ttm(i) * Freq) = (1 + r(i-1) / Freq) ^ (ttm(i-1) * Freq) * (1 + FR(i) / Freq) ^ (1)
            which is the same as 
            FR(i) = (df(i-1) / df(i) - 1) * Freq 
        """
        n = len(x)
        output = np.zeros(n)
        #output[0] = x[0][1]
        output[0] = 0.00705
        for i in range(1, n):
            output[i] = ((x[i-1][1] / x[i][1]) - 1) * Freq
        return output
        
discountFactors = bootstrap.discountCurve(100, 2, swapRates)
spotRates = bootstrap.spotRates(1, 2, np.column_stack((ttms, discountFactors)))
forwardRates = bootstrap.forwardRates(100, 2, np.column_stack((ttms, discountFactors)))

print(discountFactors)       
print("\n")       
print(spotRates)       
print("\n")       
print(forwardRates)       
            
            
 # %%
plt.plot(ttms, 100 * spotRates, label = "Spot Rate")
plt.plot(ttms, 100 * forwardRates, label = "Forward Rate")
plt.legend()
plt.ylabel("rate in %")
plt.xlabel("ttm")
# %%
