aa = [1.7 1.9 2.1 2.3];
x = 0:9;

for j = 1:length(aa)

    f = delta.^((aa(j)-3).*x);

    plot(x, f, 'LineWidth', 1.2)
    hold on
end

title('phi (1-porosity)')
xlabel('X')
ylabel('dry mass volume fraction')
legend('1.7','1.9', '2.1', '2.3')
