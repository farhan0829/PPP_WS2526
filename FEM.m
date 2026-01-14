classdef FEM < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        
    end

    methods(Static = true,...
            Access = public)
        function [K,M,Q,G,F] = buildMatricesAndVectors(fem,geo)
            %BUILDMATRCESANDVECTORS  
            % Camms OOPDE code to assembe linear system.
            [K,M,F] = fem.assema(geo,1,1,1);
            geo.makeBoundaryMatrix(geo.robinBC(0,1));
            [Q,G,~,~] = fem.assemb(geo);
              
        end

        function solveImplicit(fem,geo,setting,u0,v0,tmax,dt)
            %SOLVEIMPLICIT Implements implicite Euler method.
            
            u = u0;
            v = v0;

            [K,M,Q,G,F] = FEM.buildMatricesAndVectors(fem,geo); 
            % By s we can control the boundary conditions.  
            s = 0;
            n = geo.nPoints;
            A = [ M-dt*(setting.Buu*M-setting.Du*K+s*Q) -dt*setting.Buv*M
                   -dt*setting.Bvu*M       M-dt*(setting.Bvv*M-setting.Dv*K+s*Q)];

            C = dt*[(setting.Au*F + s*G)
                    (setting.Av*F + s*G)];

            t = 0;
            k = 0;
            while t < tmax
                % implicit Euler Step
                nl =  u.^2.*v;
                ll = [M*u
                     M*v];
                NL = dt*[setting.Cu*M*nl  
                         setting.Cv*M*nl];
                y = A\(ll+C+NL);
                u = y(1:n);
                v = y(1+n:end);
                if mod(k,setting.stepsPerFrame) == 0
                    % here a switch btw plotting u and v
                    geo.plot(u,'EdgeColor','none');                   
                    title(['t = ',num2str(t)]);                    
                    view(2);
                    drawnow
                    
                    % geo.plot(v,'EdgeColor','none');title(['t = ',num2str(t),' k = ',num2str(k)]);view(2);drawnow
                     
                end                
                t = t + dt;
                k = k + 1;
            end
        end   
        % function solveImplicit(fem,geo,coeff,u0,v0,tmax,dt)
        %     u = u0;
        %     v = v0;
        %     [K,M,Q,G,F] = FEM.buildMatricesAndVectors(fem,geo); 
        % 
        %     s = 0;
        %     A = M-dt*(-coeff.Du*K + coeff.Buu*M + coeff.Buv*M + s*Q);
        %     B = M-dt*(-coeff.Dv*K + coeff.Bvu*M + coeff.Bvv*M + s*Q);
        %     c = dt*(coeff.Au*F + s*G);
        %     d = dt*(coeff.Av*F + s*G);
        %     t = 0;
        %     k = 0;
        %     while t < tmax
        %         % implicit Euler Step
        %         y = u.^2.*v;
        %         U = A\(M*u + dt*coeff.Cu*M*y + c );
        %         V = B\(M*v + dt*coeff.Cv*M*y + d );
        % 
        %         if mod(k,5) == 0
        %             % here a switch btw plotting u and v
        %             geo.plot(U,'EdgeColor','none');
        % 
        %             title(['t = ',num2str(t),' k = ',num2str(k)]);
        % 
        %             view(2);
        %             drawnow
        % 
        %             % geo.plot(V,'EdgeColor','none');title(['t = ',num2str(t),' k = ',num2str(k)]);view(2);drawnow
        % 
        %         end
        %         u = U;
        %         v = V;
        %         t = t + dt;
        %         k = k + 1;
        %     end
        % end   
    end
end