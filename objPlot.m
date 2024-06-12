function objPlot(objList, energyList, penaltyList, totalIterations, M, K, J, fileName)
    figure;
    subplot(3, 1, 1);
    plot(1:totalIterations, objList, 'LineWidth', 2);
    title(sprintf('Objective Function Value - %s', fileName));
    xlabel('Iteration');
    ylabel('Objective Value');

    subplot(3, 1, 2);
    plot(1:totalIterations, energyList, 'LineWidth', 2);
    title('Energy Value');
    xlabel('Iteration');
    ylabel('Energy');

    subplot(3, 1, 3);
    plot(1:totalIterations, penaltyList, 'LineWidth', 2);
    title('Penalty Value');
    xlabel('Iteration');
    ylabel('Penalty');

    suptitle(sprintf('Results for M=%d, K=%d, J=%d', M, K, J));
end
