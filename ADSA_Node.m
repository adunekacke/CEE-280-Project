classdef ADSA_Node < handle

% Node class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        % 3x1 vector containing the x, y, and z coordinates of the node
        node_coord
        node_dof
        node_number
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %    Arguments
        %      node_coord:  3x1 vector containing the x, y, and z coordinates of the node
        function self = ADSA_Node(node_coord, node_number)
            self.node_coord = node_coord;
            self.node_number = node_number;
            AssignDOF(self, node_number)
        end
        
        %% Get functions of the Node Class to create copies and access the
        %copies of assigned properties from outside the class
        %Getting the coordinates of the node
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
        %Getting the DOFs of the node 
        function node_dof = GetNodeDOF(self)
            node_dof = self.node_dof;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %Assigning the DOFs to the nodes according to their node numbers
        function AssignDOF(self, node_number)
        self.node_dof=6*(node_number-1) + (1:6)';
        end
    end
end
