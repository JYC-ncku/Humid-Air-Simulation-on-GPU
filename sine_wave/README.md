Complie this code using:
```
    gcc sine_wave.c -o sine_wave.exe -lm
```
This will make a new file (Result.txt). To view in MATLAB:
```
    data = load('Result.txt')
    plot(data(:,1), data(:,2))
    xlabel('x')
    ylabel('y')
    title('y=sin(x)')
    print('Result', 'dpng', -r1000)
```

![Result.png](./Result_1000dpi.png)