# Group Assinment 1.2

*[CEGM1000 MUDE](http://mude.citg.tudelft.nl/)*

*Written by: Anna Störiko*

*Due: `<day of week>`, `<month>` `<day>`, `<year>`.*

## Part 1

**1.1 Briefly summarize how you computed the discharge from the concentration measurements and present your result.
Include the code line(s) that implement the numerical integration.**



% solution_start

We computed the discharge $Q$ from the injected solute mass $m$ and the area under the concentration curve $A$ as follows:

$$Q = \frac{m}{\int_0^{t_{\text{end}}} c(t) dt} = \frac{m}{A}$$

We approximated the integral of the concentration curve using the trapezoidal rule:

$$A = \int_0^{t_{\text{end}}} c(t) dt \approx \sum_{i=0}^{n-1} \frac{c(t_i) + c(t_{i+1})}{2} \Delta t_i$$

In the code this is implemented as follows:

```{python}
trapezoid_areas = (concentration[:-1] + concentration[1:]) / 2 * dt
total_area_trapz = np.sum(trapezoid_areas)
```

This yields an area of x mg/L·s. With an injected mass of x kg, the discharge is x L/s.

% solution_end

**1.2 <question>**

Justify your choice of numerical integration method, taking theoretical and practical considerations into account.

% solution_start
For the dataset, the trapezoidal rule is s suitable integration technique.
For a given number of integration points, it (theoretically) has a lower error compared to the left and right Riemann sum.

While Simpson's rule would have an even lower error, it requires equally spaced measurement intervals which are not available for this data set.
The error of the midpoint rule scales the same with the number of integration intervals as for the trapezoidal rule, but it requires evaluating the function at the midpoint of each interval. This is not possible with fixed measurements.
% solution_end

## Part `2, ...`

**2.1 <question>**

% solution_start
`<solution>`
% solution_end

*Copyright 2025 MUDE, TU Delft. This work is licensed under CC BY 4.0 License.*