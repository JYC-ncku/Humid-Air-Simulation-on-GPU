compile this code using:
```
gcc Perform_1D_humidity_simulation_in_Air.c -o Perform_1D_humidity_simulation_in_Air.exe -lm
```
This will make a new file (results.txt). To view in MATLAB:
```
data = load('results.txt')
plot(data(:,1), data(:,2))
xlabel('Posiition')
ylabel('H')
title('Humidity vs location')
print( 'results', 'dpng', -r1000)
```

![results.png](./results.png)