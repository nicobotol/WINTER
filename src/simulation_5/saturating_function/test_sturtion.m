x = -20:0.1:20;
y = clip(x,1,4);

figure()
plot(x, y, 'o')

function y = clip(x, a, h)
y = zeros(length(x), 1);
y = -0.5*h*(((a^2 - 2*a*x + h^2 + x.^2)/h^2).^0.5 - ((a^2 + 2*a*x + h^2 + x.^2)/h^2).^0.5);
end