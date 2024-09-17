se_mc <- function(n, m, c)
{
  se_m <- sqrt(m*(1-m)/n)
  se_c <- sqrt(c*(1-c)*((1+(n/2-1)*((1-c)/(2-c)))+((n/2-1)*c)/(1+c))/(n^2*m*(1-m)))

  c(se_m=se_m, se_c=se_c)
}



# Define the function to solve for n
solve_for_n <- function(m, c, se_m, se_c) {
  n_m <- m*(1-m)/se_m^2

  # Calculate the coefficients of the quadratic equation
  A <- se_c^2 * m * (1 - m)
  B <- -0.5 * c * (1 - c) * ((1 - c) / (2 - c) + c / (1 + c))
  C <- -c * (1 - c) * (1 - (1 - c) / (2 - c) - c / (1 + c))

  # Check if the discriminant is non-negative
  discriminant <- B^2 - 4 * A * C
  if (discriminant < 0) {
    return("No real solution")
  }

  # Use the quadratic formula to solve for n
  n1 <- (-B + sqrt(discriminant)) / (2 * A)
  n2 <- (-B - sqrt(discriminant)) / (2 * A)

  # Return the solutions
  n_c <-max(n1, n2)

  c(n_m=n_m, n_c=n_c)
}

solve_for_n(0.1, 0.75, 0.001, 0.001)

for(m in (1:50)/100)
  for(c in (51:99)/100)
  {
    n <- max(solve_for_n(m,c,0.001,0.001))
    print(paste(m,c,n))

  }


