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
    
    properties(Constant = true, Hidden = true)
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
    end
        
    methods( Static,Access = protected)
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

