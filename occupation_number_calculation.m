function occupation_number = occupation_number_calculation(MPS_tensor,mpo_tensor, M_)
% Contract the left environment, an MPS tensor, and an MPO tensor
occupation_number =1;
N= length(MPS_tensor);
 for j = 1:N    
     temp = tensorprod(occupation_number, MPS_tensor{j}, 1, 1, "NumDimensionsA", 3);
     temp = tensorprod(temp, mpo_tensor{j},[1, 3], [1, 2], "NumDimensionsA",4 );
    occupation_number = tensorprod(temp, M_{j}, [1, 3], [1, 2], "NumDimensionsA",4);
 end
