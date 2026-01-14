%% lagrange12D
% Class that implements Lagrange-1 elements in 2D.

%%Copy style%
%handle

%%Inheritance
% lagrange22D < finiteElements2D < finiteElements < handle


classdef Lagrange12D < finiteElements2D
    % class for build Lagrange L1 FE in 2D
    % Inherits everything from abstract class finiteElements
    % Only different Matrices were defined as constant props.
    
    
    properties(Constant = true, Access = public)
        %% Constant properties with Access = public
        %
        % * S1, S2 S3   (double) *vectors* that store the element stiffness matrices.
        % * C1 C2       (double) *vectors* that store the element gradient matrices
        % * M           (double) *vector* that stores the element mass matrix
        % * F           (double) vector that stores the source mass vector.
        % * idx         (double) vector that stores order of Points in element field.
        % * boundaryIndex (double) = [1 2]
        S1 = 0.5*[ 1
                  -1
                   0
                  -1
                   1
                   0
                   0
                   0
                   0];
               
        S2 = 0.5*[ 2 
                  -1 
	              -1 
                  -1 
                   0 
                   1  
                  -1
                   1
                   0];
              
        S3 = 0.5*[ 1  
                   0  
	              -1  
                   0 
                   0  
                   0  
                  -1
                   0
                   1];
              
        % mass
        M = 1/24*[ 2  
	               1  
	               1  
                   1   
                   2  
                   1  
                   1
                   1 
                   2];
               
        % RHS vector
        F = 1/6*[ 1
                  1
                  1]; 
              
        % Matrices convection terms 
        % for use by "sparse" in vector form
        C1 = [-1/6   
           -1/6   
           -1/6  
            1/6
            1/6
            1/6
            0
            0
            0];
        
        C2 = [-1/6   
           -1/6   
           -1/6 
            0
            0
            0
            1/6
            1/6
            1/6]; 
        
        % index vector for selecting coordinates
        idx  = 1:3;          
        boundaryIndex = [1,2];        
    end
    
    properties(Constant = true, Access = protected)
        %% Constant properties wiht Access = protected
        %
        % boundaryElements@finiteElements = lagrange11D
        %
        
        boundaryElements = Lagrange11D
    end
    
    methods(Static = true, Access=public)  
        %% Overwritten method
        %
        % * makeIndex
        
        % 
        % This method computes the explicite given structure of 
        % point-element relation. We need it in the abstract methods
        % iherited from superclasses.
        %
        function [idxvec0,idxvec1,idxvec2] = makeIndex(idx,np)
            idxvec0 = reshape(idx,1,np*3);
            idxvec1 = reshape([idx;idx;idx],1,np*9);
            idxvec2 = reshape([idx(1,:);idx(1,:);idx(1,:);...
                                idx(2,:); idx(2,:); idx(2,:);...
                                idx(3,:);idx(3,:);idx(3,:)],1,9*np);     
        end  
    end 
    
    methods(Access=public)
        %% Methods with Access = public
        %
        function [DX,DY] = gradientMatrices(obj,grid)
            %%
            % * gradientMatrices IN: gridd OUT: double,double
            % 
            % Method that computes matrices DX DY, such that 
            %
            % grad u = [(DX*u)' ,(DY*u)']                       
            %
            % at the center of all triangles. 
            % DX and DY are nElements x nPoints matrices. Note that this is
            % exact for linear FE's but an approximation with only
            % linear convergence for the function u. 
            
            % We want to use makeJ, so it is not static...
            
            % The indices of the first, second and third points in each
            % triangle. 
            idx1 = grid.t(1,:);
            idx2 = grid.t(2,:);
            idx3 = grid.t(3,:);
            
            % p1 are the  values of the "first" point, p2 the values of
            % the "second" point p3 the values of the third point 
            % of all triangles. (2 x nElement matrices)
            p1 = grid.p(:,(idx1));
            p2 = grid.p(:,(idx2));
            p3 = grid.p(:,(idx3));
            
            % Compute Jacobi determinat  
            J = obj.makeJ(grid);
            
            
            p12 = (p1(2,:)-p2(2,:))./J;
            p31 = (p3(2,:)-p1(2,:))./J;
            p23 = (p2(2,:)-p3(2,:))./J;
           
            DX = sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3],[p23 p31 p12] ,...
                grid.nElements,grid.nPoints);
            
            p21 = (p2(1,:)-p1(1,:))./J;
            p13 = (p1(1,:)-p3(1,:))./J;
            p32 = (p3(1,:)-p2(1,:))./J;
            
            DY = sparse([1:grid.nElements,1:grid.nElements,1:grid.nElements],...
                [idx1 idx2 idx3],[p32 p13 p21] ,...
                grid.nElements,grid.nPoints);
        end  
    end
end
