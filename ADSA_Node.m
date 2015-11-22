classdef ADSA_Node < handle

%Node class for a three-dimentional elastic frame structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Functions Called
%     Public:
%       ADSA_Node      Constructs the node object and stores its
%                        properties
%       GetNodeCoord   Getter function to get the nodal coordinates of any
%                        node
%       GetNodeDOF     Getter function to get the nodal Degrees of Freedom
%                        of any node
%     Private:
%       AssignDOF      Assign the DOFs of a node according to its node
%                        number and saves it as a property
%       
%
%  Dictionary of Variables
%     Input Information and Properties:
%       nodeCoord      ==  3x1 vector containing the x, y, and z 
%                            coordinates of the node            
%       nodeDOF        ==  6x1 vector containing the 6 Degrees of Freedom 
%                            (displacement along x, y and z directions and 
%                            rotation about x, y and z directions
%                            respectively) of the node corresponding to its
%                            node number
%       nodeNumber     ==  The number assigned to the node which is used as
%                            its identifier and to assign its DOFs
%
%     Local Variables:
%       No local variables used in the Node Class
    
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
