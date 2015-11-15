classdef ADSA_Node < handle

% Node class for a 3-dimensional framed structure
    
    %% Private properties go here
    properties (Access = private)
        % 3x1 vector containing the x, y, and z coordinates of the node
        nodeCoord
        nodeDOF
        nodeNumber
    end
    
    %% Public methods go here
    methods (Access = public)
        
        %% Constructor
        function self = ADSA_Node(nodeCoord, nodeNumber)
            self.nodeCoord = nodeCoord;
            self.nodeNumber = nodeNumber;
            AssignDOF(self)
        end
        
        
        %% Get functions of the Node Class to create copies and access the
        %copies of assigned properties from outside the class
        
        %Getting the coordinates of the node
        function nodeCoord = GetNodeCoord(self)
            nodeCoord = self.nodeCoord;
        end
        
        
        %Getting the DOFs of the node 
        function nodeDOF = GetNodeDOF(self)
            nodeDOF = self.nodeDOF;
        end
        
    end
    
    
    %% Private methods go here
    methods (Access = private)
        
        %Assigning the DOFs to the nodes according to their node numbers
        function AssignDOF(self)
            
        self.nodeDOF=6*(self.nodeNumber-1) + (1:6)';
        
        end
        
    end
    
end
