# Generalized Penalty Function Approach

We used a new identification strategy to identify the aggregate demand and the aggregate supply shock on real U.S. data. The simulation study and the robust test show that this identification strategy works well

* the estimated implied impulse response functions are consistent with the model-implied impulse response functions
* the estimated shocks and actuals shocks are very highly correlated. 

The empirical results indicate that GDP is less persistent than inflation in response to the demand shock, whereas GDP is more persistent than inflation to the supply shock.



### Identification Scheme

In order to identify the demand and supply shocks, according to the AD-AS model, we have the following information at hand for identification: the demand shock, $\epsilon^D_t$, imposes a positive impact on both output, $Y_t$, and inflation, $\pi_t$, and the supply shock, $\epsilon^S_t$, imposes a positive impact on the output, $Y_t$, and a negative impact on the inflation, $\pi_t$.

The response of output (inflation) to demand shocks should be increasing functions of the size of the demand shock and the response of output (inflation) to the supply shock $\epsilon^S_t$ should be an increasing (a decreasing)  function of the size of the supply shock.




### Impulse responses

We use a simple DSGE model to simulate data and use the simulated data to estimate the impulse response.

![baselineCompD](graph/baselineCompD.png)

![baselineCompS](graph/baselineCompS.png)



### What's more...

We also extend our method to bayesian estimation, codes can be found in the folder: bayesian.

<img src="graph/demand shock.png" alt="demand shock" style="zoom:80%;" />

<img src="./graph/supply shock.png" alt="supply shock" style="zoom:80%;" />



### Thanks

This is a joint paper with my advisor, Dr. Marco Brianti. I wish to extend my deepest thanks to him for providing generous guidance on this project and the SVAR model's identification and estimation techniques.
