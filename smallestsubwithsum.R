smallestSubWithSum = function(arr, n, x){  
  
  # Initialize current sum and minimum length 
  curr_sum = 0 
  min_len = n + 1 
  
  # Initialize starting and ending indexes 
  start = 1
  solution_start = 1
  end = 1
  solution_end = 1
  
  while(end < n){ 
    # Keep adding array elements while current sum is smaller than x 
    while(curr_sum <= x & end < n){ 
      #print("end")
      # Ignore subarrays with negative sum if x is positive. 
      if(curr_sum <= 0 & x > 0){ 
        start = end 
        curr_sum = 0 
      }
      curr_sum = curr_sum + arr[end] 
      end = end + 1
      
    }
    
    # If current sum becomes greater than x. 
    while(curr_sum > x & start < n){ 
      # Update minimum length if needed 
      if(end - start < min_len){  
        min_len = end - start
        solution_start = start
        solution_end = end-1
      }
      # remove starting elements 
      curr_sum = curr_sum - arr[start] 
      start = start + 1
    }
    
  }
  return(c(solution_start, solution_end))
  }

