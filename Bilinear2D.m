classdef Bilinear2D < finiteElements2D
    % class definition of bilinear elements 
    % on rectangular meshes. Note this definition works only for
    % elements with parallel edges.
    
    properties(Constant = true, Access = public)
        S1 = 1/6 * reshape([ 2	-2	-1	 1	
                            -2	 2	 1	-1	
                            -1	 1   2	-2	
                             1	-1	-2	 2],16,1) 
             
        S2 = 1/2 * reshape([ 1   0  -1   0
                             0  -1   0   1
                            -1   0   1   0
                             0   1   0  -1],16,1)  
        
        S3 = 1/6 * reshape([ 2   1  -1  -2
                             1   2  -2  -1
                            -1  -2   2   1
                            -2  -1   1   2],16,1)  
        
        M = 1/36 * reshape([ 4	 2	 1	 2	
                             2   4   2   1
                             1   2   4   2 
                             2   1   2   4],16,1)  
                         
        F = 1/4 * [ 1
                    1
                    1
                    1]
    end
    
    properties(Constant = true, Access = public)
        % The order of point indeces in the e field. Trivialy [1,2] for P1 
        % elements but importand for extended meshes.
        boundaryIndex = 1:2
    end
    
    properties(Constant = true, Access = protected)
        boundaryElements = Lagrange11D
    end
    
    properties(Constant = true,...
            Hidden = true)
        % Non trivial order, see the book of Schwarz. 
        idx = [1 2 4 3] 
    end
    
    methods(Access = public, Static = true)
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            idxvec0 = reshape(idx,1,np*4);
            idxvec1 = reshape([idx;idx;idx;idx],1,np*16);
            idxvec2 = reshape([idx(1,:);idx(1,:);idx(1,:);idx(1,:);...
                            idx(2,:); idx(2,:); idx(2,:); idx(2,:);...
                            idx(3,:);idx(3,:); idx(3,:);idx(3,:);...
                            idx(4,:);idx(4,:); idx(4,:);idx(4,:)],1,16*np);     
        end  




        function H = periodic(Q)
            n1 = Q.e(2,Q.e(end,:) == 1);
            % n2 = Q.e(2,Q.e(end,:) == 2);
            % n3 = Q.e(2,Q.e(end,:) == 3);
            n4 = Q.e(2,Q.e(end,:) == 4);
            m = nnz(Q.e(end,:) == 2)+1;
            n = nnz(Q.e(end,:) == 1)+1;
             
            H = spalloc(Q.nBoundaryPoints,m*n,Q.nBoundaryPoints*4);
            
            k1 = 1;
            
            % left lower corner
            H(k1,1) = 1;
            H(k1,m*n) = -1;
            k1 = k1 + 1;
            
            for k2 = 1:length(n1)-1
                H(k1,n1(k2))=1;
                H(k1,n1(k2)+m-1)= -1; 
                k1 = k1 + 1;
            end
            for k2 = length(n4)-1:-1:1
                H(k1,n4(k2))=1;
                H(k1,n4(k2)+m*(n-1))= -1; 
                k1 = k1 + 1;
            end
            % corner points diagonal symmetric!
            
            H(k1,m) = 1;
            H(k1,m*(n-1)+1) = -1;
            k1 = k1 + 1;
            hy = mean(Q.ygrid.h);
            % % first order derivative
            for k2 = 1:length(n1)-1
                H(k1,n1(k2))=-1/hy;
                H(k1,n1(k2)+1)= 1/hy;
                H(k1,n1(k2)+m-2)= -1/hy; 
                H(k1,n1(k2)+m-1)= 1/hy; 
                k1 = k1 + 1;
            end
            hx = mean(Q.xgrid.h);
            for k2 = length(n4)-1:-1:1
                H(k1,n4(k2))=-1/hx;   
                H(k1,n4(k2)+m)= 1/hx; 
            
                H(k1,n4(k2)+m*(n-2))= 1/hx;
                H(k1,n4(k2)+m*(n-1)) = -1/hx; 
                k1 = k1 + 1;
            end
            h = sqrt(hx^2+hy^2);
            % Derivatives in corners
            H(k1,1) = -1/h;
            H(k1,m+2) = 1/h; 
            H(k1, m*n) = 1/h;
            H(k1, m*(n-1)-1) = -1/h;
             
            k1 = k1 + 1;
            
            H(k1,m) = -1/h;
            H(k1,2*m-1) = 1/h; 
            H(k1, (n-1)*m+2) = -1/h;
            H(k1, m*(n-2)+1) = 1/h;
            
        end

    end
        
    methods( Static = true,...
            Access = protected)
        % !!!! Dummies !!!! If you don't want to implement adaptive error
        % control, you can leave it as it is. Otherwise implements this
        % functions.
        function fluxThrougElementEdges= fluxThroughEdges(gridObj,u,c)
        end
        function normFminusau = localErrorL2(gridObj,a,f) 
        end
        function jumps = fluxJumps(gridObj,fluxThroughElementEdges,order)
        end
    end
    
end

