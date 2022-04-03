%Forward Backward
%delta is the distance between two consecutive range measurement 
function p=forback1(D,delta,VLmax,TI)
dmax=2/3*VLmax*TI';
deltanum=[1 2 4 8]*double(delta>dmax)+1;

transition=zeros(2,2,16);
transition(:,:,1)=[0.5 0.5;0.5 0.5];  %delta(1,i)<dmax && delta(2,i)<dmax && delta(3,i)<dmax && delta(4,i)<dmax (1)
transition(:,:,2)=[0.1 0.9;0.5 0.5];      %delta(1,i)>dmax && delta(2,i)<dmax && delta(3,i)<dmax && delta(4,i)<dmax (2)
transition(:,:,3)=[0.9 0.1;0.5 0.5];      %delta(1,i)<dmax && delta(2,i)>dmax && delta(3,i)<dmax && delta(4,i)<dmax (3)
transition(:,:,4)=[0.5 0.5;0.5 0.5];      %delta(1,i)>dmax && delta(2,i)>dmax && delta(3,i)<dmax && delta(4,i)<dmax (4)
transition(:,:,5)=[0.5 0.5;0.1 0.9];      %delta(1,i)<dmax && delta(2,i)<dmax && delta(3,i)>dmax && delta(4,i)<dmax (5)
transition(:,:,6)=[0.1 0.9;0.1 0.9];          %delta(1,i)>dmax && delta(2,i)<dmax && delta(3,i)>dmax && delta(4,i)<dmax (6)
transition(:,:,7)=[0.9 0.1;0.1 0.9];          %delta(1,i)<dmax && delta(2,i)>dmax && delta(3,i)>dmax && delta(4,i)<dmax (7)
transition(:,:,8)=[0.5 0.5;0.1 0.9];          %delta(1,i)>dmax && delta(2,i)>dmax && delta(3,i)>dmax && delta(4,i)<dmax (8)
transition(:,:,9)=[0.5 0.5;0.9 0.1];      %delta(1,i)<dmax && delta(2,i)<dmax && delta(3,i)<dmax && delta(4,i)>dmax (9)
transition(:,:,10)=[0.1 0.9;0.9 0.1];         %delta(1,i)>dmax && delta(2,i)<dmax && delta(3,i)<dmax && delta(4,i)>dmax(10)
transition(:,:,11)=[0.9 0.1;0.9 0.1];         %delta(1,i)<dmax && delta(2,i)>dmax && delta(3,i)<dmax && delta(4,i)>dmax (11)
transition(:,:,12)=[0 0;0.9 0.1];         %delta(1,i)>dmax && delta(2,i)>dmax && delta(3,i)<dmax && delta(4,i)>dmax (12)
transition(:,:,13)=[0.5 0.5;0.5 0.5];     %delta(1,i)<dmax && delta(2,i)<dmax && delta(3,i)>dmax && delta(4,i)>dmax (13)
transition(:,:,14)=[0.1 0.9;0.5 0.5];         %delta(1,i)>dmax && delta(2,i)<dmax && delta(3,i)>dmax && delta(4,i)>dmax (14)
transition(:,:,15)=[0.9 0.1;0.5 0.5];         %delta(1,i)<dmax && delta(2,i)>dmax && delta(3,i)>dmax && delta(4,i)>dmax (15)
transition(:,:,16)=[0.5 0.5;0.5 0.5]; %delta(1,i)>dmax && delta(2,i)>dmax && delta(3,i)>dmax && delta(4,i)>dmax (16)

observation=zeros(2,2,16);
for i=1:16
    observation(:,:,i)=[0.5 0;0 0.5];
end
    
if D==1
   observation(:,:,2)=[0.8 0;0 0.2];
elseif D==2
   observation(:,:,2)=[0.2 0;0 0.8];
end        


f=zeros(2,size(delta,2)+1);
  for i=1:size(delta,2)
      trans=transition(:,:,deltanum(i));
      obs=observation(:,:,deltanum(i));
   if i==1 
       f(:,1)=[obs(1,1);obs(2,2)];
   end
   if f(1,i)>=f(2,i) && obs(1,1)==0.5 %if Observation is 50%, assume state did not change
       obs=obs*[1.3 0; 0 1/1.3];
   elseif f(1,i)<f(2,i) && obs(1,1)==0.5
       obs=obs*[1/1.3 0; 0 1.3];
   end
   f(:,i+1)=obs*trans*f(:,i)/sum(obs*trans*f(:,i));
  end
  
  b=zeros(2,size(delta,2)+1);
  for i=size(delta,2):-1:1
      trans=transition(:,:,deltanum(i));
      obs=observation(:,:,deltanum(i));
   if i==size(delta,2) 
       b(:,size(delta,2)+1)=[obs(1,1);obs(2,2)];
   end
   if b(1,i+1)>=b(2,i+1) && obs(1,1)==0.5 %if Observation is 50%, assume state did not change
       obs=obs*[1.3 0; 0 1/1.3];
   elseif b(1,i+1)<b(2,i+1) && obs(1,1)==0.5
       obs=obs*[1/1.3 0; 0 1.3];
   end
   b(:,i)=trans*obs*b(:,i+1)/sum(obs*trans*b(:,i+1));
  end
  
 p=f.*b./[sum(f.*b);sum(f.*b)];  