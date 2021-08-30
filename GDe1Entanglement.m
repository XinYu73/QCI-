function  outputArg1  = GDe1Entanglement(p,Theta,Phi,y,measurement)
format long
Header=1;
        if p==0&&measurement=="Negativity"
            outputArg1=2.*(y.*sqrt(3+cos(Theta))).*(sin(Theta./2).^2)./(sqrt(2).*(1+y.^2+2.*y.*cos(Theta./2).^3.*cos(Phi)));%for Gamma,for array of coefs
            return ;
        end
%--------------------------------
        i=sqrt(-1);
        rho11=(-16*(-2+p)^3+16*y*(-y*(-1+(-1+p)*cos(Theta))^3-2*(-2+p)^3*cos(Theta/2)^3*cos(Phi)))/(128*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho12=-((exp(-i*Phi)*(-1+p)*y*((-2+p)^2*cos(Theta/2)+exp(i*Phi)*y*(1+cos(Theta)-p*cos(Theta))^2)*sin(Theta))/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho14=(exp(-i*Phi)*(-1+p)^2*y*cos(Theta/2)*(2-p-exp(i*Phi)*y*cos(Theta/2)*(-1+(-1+p)*cos(Theta)))*sin(Theta/2)^2)/(2*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho15=-(((-1+p)*y*(exp(i*Phi)*(-2+p)^2*cos(Theta/2)+y*(1+cos(Theta)-p*cos(Theta))^2)*sin(Theta))/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho16=-(((-1+p)^2*y^2*(-1+(-1+p)*cos(Theta))*sin(Theta)^2)/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho18=-(((-1+p)^3*y^2*sin(Theta)^3)/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho22=(4*(-2+p)^2*p-2*(-1+(-2+p)*p)*y^2+y*((1+p*(5+3*(-3+p)*p))*y*cos(Theta)+(-1+p)^2*y*(-2*cos(2*Theta)+(-1+p)*cos(3*Theta))+8*(-2+p)^2*p*cos(Theta/2)^3*cos(Phi)))/(32*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho24=(exp(-i*Phi)*(-1+p)*y*(2*(-2+p)*p*cos(Theta/2)+exp(i*Phi)*y*(-1+(-2+p)*p+(-1+p)^2*cos(2*Theta)))*sin(Theta))/(16*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho25=((-1+p)^2*y*cos(Theta/2)*(-exp(i*Phi)*(-2+p)+y*cos(Theta/2)*(1+cos(Theta)-p*cos(Theta)))*sin(Theta/2)^2)/(2*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho26=((-1+p)*y*(2*exp(i*Phi)*(-2+p)*p*cos(Theta/2)+y*(-1+(-2+p)*p+(-1+p)^2*cos(2*Theta)))*sin(Theta))/(16*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho28=((-1+p)^2*y^2*(1+(-1+p)*cos(Theta))*sin(Theta)^2)/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho44=(2*y^2-2*(-2+p)*p*(2*p+y^2)+y*((-1+p)*y*(cos(Theta)-3*(-2+p)*p*cos(Theta)-2*(-1+p)*cos(2*Theta)-(-1+p)^2*cos(3*Theta))-8*(-2+p)*p^2*cos(Theta/2)^3*cos(Phi)))/(32*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho45=-(((-1+p)^3*y*(exp(i*Phi)+y*cos(Theta/2)^3)*sin(Theta/2)^3)/(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho46=((-1+p)^2*y*(4*exp(i*Phi)*p*cos(Theta/2)+y*(1+p+2*p*cos(Theta)+(-1+p)*cos(2*Theta)))*sin(Theta/2)^2)/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho48=-(((-1+p)*y*(exp(i*Phi)*p^2*cos(Theta/2)+y*(1+(-1+p)*cos(Theta))^2)*sin(Theta))/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho54=-((exp(-i*Phi)*(-1+p)^3*y*(1+exp(i*Phi)*y*cos(Theta/2)^3)*sin(Theta/2)^3)/(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho58=(exp(-i*Phi)*(-1+p)^2*y*cos(Theta/2)*(p+exp(i*Phi)*y*cos(Theta/2)*(1+(-1+p)*cos(Theta)))*sin(Theta/2)^2)/(2*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        rho68=-((exp(-i*Phi)*(-1+p)*y*(p^2*cos(Theta/2)+exp(i*Phi)*y*(1+(-1+p)*cos(Theta))^2)*sin(Theta))/(8*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi))));
        rho88=(4*p^3+2*(5+3*(-2+p)*p)*y^2+y*((-1+p)*y*(3*(5+(-2+p)*p)*cos(Theta)+(-1+p)*(6*cos(2*Theta)+(-1+p)*cos(3*Theta)))+8*p^3*cos(Theta/2)^3*cos(Phi)))/(32*(1+y^2+2*y*cos(Theta/2)^3*cos(Phi)));
        %---------------------------------ppt(epsilon(rho0))
        rho=[rho11	rho12	rho12	rho14	rho15	rho16	rho16	rho18;
            rho15	rho22	rho16	rho24	rho25	rho26	rho18	rho28;
            rho15	rho16	rho22	rho24	rho25	rho18	rho26	rho28;
            rho25	rho26	rho26	rho44	rho45	rho46	rho46	rho48;
            rho12	rho14	rho14	rho54	rho22	rho24	rho24	rho58;
            rho16	rho24	rho18	rho58	rho26	rho44	rho28	rho68;
            rho16	rho18	rho24	rho58	rho26	rho28	rho44	rho68;
            rho18	rho28	rho28	rho68	rho46	rho48	rho48	rho88;];
%---------------------------------------------------
switch measurement
    case "Negativity"%--------------------------------
        if p==0
            outputArg1=2.*(y.*sqrt(3+cos(Theta))).*(sin(Theta./2).^2)./(sqrt(2).*(1+y.^2+2.*y.*cos(Theta./2).^3.*cos(Phi)));%for Gamma,for array of coefs
            return ;
        elseif p==Header
            outputArg1= -1;
            return ;
        end
        eigValues=real(eig(rho));% possiblely wrong
        temp=0;
        for j=1:length(eigValues)
            if eigValues(j)<0
                temp=temp+eigValues(j);
            end
        end
        if temp<0
            outputArg1=-2*temp;
            return ;
        else
            outputArg1=-1;
            return ;
        end
    case "Concurrence"%--------------------------------
        if p==Header
            outputArg1= -1;
            return ;
        end
        rhoBC=rho(1:size(rho,1)/2,1:size(rho,1)/2)+rho(size(rho,1)/2+1:size(rho,1),size(rho,1)/2+1:size(rho,1));
        sigmay=[0,-i;i,0];
        rhocon=rhoBC*kron(sigmay,sigmay)*conj(rhoBC)*kron(sigmay,sigmay);
        lambda=sort(real(sqrt(eig(rhocon))),'descend');
        if p==0
            outputArg1=2*lambda(1)-sum(lambda);
            return;
        end
        if 2*lambda(1)-sum(lambda)>0
            outputArg1=2*lambda(1)-sum(lambda);
        else
            outputArg1=-1;
        end
    case "Witness"%--------------------------------
        if p==Header
            outputArg1= -1;
            return ;
        end
        rho0=inversePaitialTranspose(rho);
        rho0=CB2blocks(rho0,3);
        outputArg1=-pptmixerPI(rho0,3);%to cope with other measurement
    case "ThreeTangle"%--------------------------------
        if p==Header
            outputArg1= -1;
            return ;
        end
        rho0=inversePaitialTranspose(rho);
        t_res=NaN;
        while isnan(t_res)
        [chi, lambda] = densityEig(rho0);
        % Convex-sum function handles for the tangle
        f_cr = @(x) convexSum(x, @tangle, ...
                        chi, lambda);
        g_cr = @(x) grad_convexSum(x, @tangle, ...
                        @grad_tangle, chi, lambda);   
        % Dimensions of the Stiefel manifold
        r = rank(rho0);
        k = r+4 ;
        % Objective function and gradient
        f_opt = @(x) f_cr(buildUnitary_script(x, k, r));
        g_opt = @(x) grad_eh_adapt(x, k, r, g_cr);
%       [f_eh, g_eh] = createEHFunctions(rho0, k, r, @tangle, @grad_tangle);
        X0 = 2*pi*randn(1, dimSt(k, r));
        [t_res, ~, ~] = bfgs_min(f_opt, g_opt, X0);
        if isnan(t_res)
            disp("NaN caught")
        end
        end
        outputArg1=t_res;
end
end
function rho0=inversePaitialTranspose(rhoppt)
rho0temp=zeros(8,8);
rho0temp(1:4,1:4)=rhoppt(1:4,1:4);
rho0temp(5:8,5:8)=rhoppt(5:8,5:8);
rho0temp(5:8,1:4)=rhoppt(1:4,5:8);
rho0temp(1:4,5:8)=rhoppt(5:8,1:4);
rho0=rho0temp;
end