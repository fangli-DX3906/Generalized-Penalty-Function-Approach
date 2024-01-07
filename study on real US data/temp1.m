clc
clear
radius=1;
simNum=25;
pi_est=[];


for i=100:100:10000
    pi_sum=0;
    for j=1:simNum
        dots=rands(i,2);
        distance=dots(:,1).^2+dots(:,2).^2;
        check=sqrt(distance)<=radius;
        pi_sum=pi_sum+4*sum(check)/i;
    end
    pi_est=[pi_est pi_sum/simNum];
end

figure
plot(1:size(pi_est,2), pi_est)
hold on
plot(1:size(pi_est,2), pi*ones(1,size(pi_est,2)),'b-')
ylim([3,3.2])