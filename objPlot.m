function objPlot(objList, energyList, penaltyList, totalIterations, M, K, J, fileName)
    figure;
    subplot(3, 1, 1);
    plot(1:totalIterations, objList, 'LineWidth', 2);
    title(sprintf('Objective Function Value - %s', fileName),'FontSize', 24);
    xlabel('Iteration','FontSize', 24);
    ylabel('Objective Value','FontSize', 24);

    subplot(3, 1, 2);
    plot(1:totalIterations, energyList, 'LineWidth', 2);
    title('Energy Value','FontSize', 24);
    xlabel('Iteration','FontSize', 24);
    ylabel('Energy','FontSize', 24);

    subplot(3, 1, 3);
    plot(1:totalIterations, penaltyList, 'LineWidth', 2);
    title('Penalty Value','FontSize', 24);
    xlabel('Iteration','FontSize', 24);
    ylabel('Penalty','FontSize', 24);

    suptitle(sprintf('Results for M=%d, K=%d, J=%d', M, K, J));
end
