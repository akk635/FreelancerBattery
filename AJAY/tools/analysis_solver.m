
D = size(Ut);
timestep = 1:D(3);

figure(11)

subplot(2,1,1)
plot(timestep,P.residual_U,'-b','LineWidth',1)
xlabel('TIME STEP')
ylabel('RESIDUAL')

subplot(2,1,2)
plot(timestep,P.iteration_U,'-r','LineWidth',1)
xlabel('TIME STEP')
ylabel('ITERATION No.')