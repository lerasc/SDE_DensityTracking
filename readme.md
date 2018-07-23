The implementation consists of one single class that tracks the density distribution of the solution to a stochastic differential equation (SDE). 
Starting from a generic SDE
<html>
<head>
<title>LaTeX4Web 1.4 OUTPUT</title>
<style type="text/css">
<!--
 body {color: black;  background:"#FFCC99";  }
 div.p { margin-top: 7pt;}
 td div.comp { margin-top: -0.6ex; margin-bottom: -1ex;}
 td div.comb { margin-top: -0.6ex; margin-bottom: -.6ex;}
 td div.norm {line-height:normal;}
 td div.hrcomp { line-height: 0.9; margin-top: -0.8ex; margin-bottom: -1ex;}
 td.sqrt {border-top:2 solid black;
          border-left:2 solid black;
          border-bottom:none;
          border-right:none;}
 table.sqrt {border-top:2 solid black;
             border-left:2 solid black;
             border-bottom:none;
             border-right:none;}
-->
</style>
</head>
<body>

<table cellspacing=0  border=0 align=center>
<tr>
  <td nowrap align="center">
    
	\dd X<sub>t</sub> = <font face=symbol>m</font>(X<sub>t</sub>)&nbsp; \dd t + <font face=symbol>s</font>(X<sub>t</sub>) &nbsp;\dd W<sub>t</sub>
<a name="eq0">&nbsp;&nbsp;&nbsp;&nbsp;<font color=blue>(0)</font>
  </td>
</tr>
</table>
</body>
</html>
with some time-independent drift and volatility functions $\mu$ and $\sigma$, and initial position $X_0$, 
the class calculates the probability density $p(x,t)$, to be at position $x$ at time $t$. 
The class can furthermore deal with absorbing or reflective boundary conditions.

See description.pdf and documentation inside SDE_DensityTracking.py for details. 