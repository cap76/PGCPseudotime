function K = covBranchingProcess_v5c(hyp, x, z, i)

%A two component branching process.

if nargin<2, K = '10'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

k1={'covMaterniso',3}; %Base kernel
 
location1  = hyp(1);
steepness1 = hyp(2);
location2  = hyp(3);
steepness2 = hyp(4);

covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
%covN2 = {@covSum,{covSEiso,{@covChangePointMultiD, {1, @covZero, @covSEiso}}}};
%covN3 = {@covSum,{covSEiso,{@covSum,{{@covChangePointMultiD, {1, @covZero, @covSEiso}},{@covChangePointMultiD, {1, @covZero, @covSEiso}}}}}};

if nargin<4 %Covariances  
    ind1 = find(x(:,2)==1);
    ind2 = find(x(:,2)==2);
    
    if dg    
        %K = zeros(size(x,1),size(x,1));
        K   = feval(k1{:},[hyp(9),hyp(10)],x(:,1));
        K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1));
        K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1));
        %K(ind1,ind1) = 
        K(ind1,ind1) = K(ind1,ind1)+K2;
        K(ind2,ind2) = K(ind2,ind2)+K3;
        %K    = [K1(ind1,ind1),K1(ind1,ind2);
        %        K1(ind2,ind1),K1(ind2,ind2)+K2];        
        K    = diag(K);
    else        
    if  xeqz   
    K    = feval(k1{:},[hyp(9),hyp(10)],x(:,1));
    %covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
    K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1));
    K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1));    
    %K = [K1(ind1,ind1),K1(ind1,ind2);
    %     K1(ind2,ind1),K1(ind2,ind2)+K2(ind2,ind2)];
    K(ind1,ind1) = K(ind1,ind1)+K2;
    K(ind2,ind2) = K(ind2,ind2)+K3;
    else %Cross covariances

        ind4 = find(z(:,2)==1);
        ind5 = find(z(:,2)==2);

         
        K    = feval(k1{:},[hyp(9),hyp(10)],x(:,1),z(:,1));
        %covN1 = {@covChangePointMultiD, {1, @covZero, {'covMaterniso',3}}};
        %covN2 = {@covChangePointMultiD, {1, @covZero, @covSEiso}};
        %K = zeros(size(x,1),1);
        K2   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),z(ind4,1));
        K3   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),z(ind5,1));        

    try
    K(ind1,ind4) = K(ind1,ind4)+K2;
    end        
        
    try
    K(ind2,ind5) = K(ind2,ind5)+K3;
    end
         %K = [K1(ind1,ind4),K1(ind1,ind5);            
         %     K1(ind2,ind4),K1(ind2,ind5)+K2(ind2,ind5)];
        
    end    
end
    
else
    
ind1 = find(x(:,2)==1);
ind2 = find(x(:,2)==2);
  % derivatives
  if i==1
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),1);      
  elseif i==2
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),2);          
  elseif i==3
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),1);      
  elseif i==4
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),2);          
  elseif i==5
    K = zeros(size(x,1),size(x,1));
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),3);
  elseif i==6
    K = zeros(size(x,1),size(x,1));      
    K(ind1,ind1)   = feval(covN1{:}, [location1,steepness1,hyp(5),hyp(6)], x(ind1,1),x(ind1,1),4);    
  elseif i==7
    K = zeros(size(x,1),size(x,1));
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),3);
  elseif i==8
    K = zeros(size(x,1),size(x,1));      
    K(ind2,ind2)   = feval(covN1{:}, [location2,steepness2,hyp(7),hyp(8)], x(ind2,1),x(ind2,1),4);    
  elseif i==9
     K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),x(:,1),1);
  elseif i==10
     K = feval(k1{:},[hyp(9),hyp(10)],x(:,1),x(:,1),2);
  else
    error('Unknown hyperparameter')
  end
end