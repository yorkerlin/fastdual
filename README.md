<B><FONT COLOR="#663300">%</FONT></B>example one
clear all<B><FONT COLOR="#663300">;</FONT></B> close all<B><FONT COLOR="#663300">;</FONT></B>

x1<B><FONT COLOR="#663300">=[</FONT></B><FONT COLOR="#996600">0.8822936</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#996600">0.7160792
0.9178174</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#996600">0.0135544</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#996600">0.5275911</FONT><B><FONT COLOR="#663300">];</FONT></B>

x2<B><FONT COLOR="#663300">=[-</FONT></B><FONT COLOR="#996600">0.9597321
0.0231289
0.8284935
0.0023812</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#996600">0.7218931</FONT><B><FONT COLOR="#663300">];</FONT></B>

x<B><FONT COLOR="#663300">=[</FONT></B>x1' x2'<B><FONT COLOR="#663300">];</FONT></B>

y<B><FONT COLOR="#663300">=[</FONT></B><FONT COLOR="#999900">1</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#999900">1
1</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#999900">1</FONT><B><FONT COLOR="#663300">
-</FONT></B><FONT COLOR="#999900">1</FONT><B><FONT COLOR="#663300">];</FONT></B>

hyp<B><FONT COLOR="#663300">.</FONT></B>cov<B><FONT COLOR="#663300"> =</FONT></B> log<B><FONT COLOR="#663300">([</FONT></B><FONT COLOR="#999900">2</FONT><B><FONT COLOR="#663300">,</FONT></B><FONT COLOR="#999900"> 2</FONT><B><FONT COLOR="#663300">]);</FONT></B>
hyp<B><FONT COLOR="#663300">.</FONT></B>lik<B><FONT COLOR="#663300">=[];</FONT></B>
cov<B><FONT COLOR="#663300"> = {</FONT></B><FONT COLOR="#009900">'covSEiso'</FONT><B><FONT COLOR="#663300">};</FONT></B>
lik<B><FONT COLOR="#663300"> = {</FONT></B>@likLogistic<B><FONT COLOR="#663300">}

[</FONT></B>a b c d e<B><FONT COLOR="#663300">]=</FONT></B>approxDualWithLBFGSForLogit<B><FONT COLOR="#663300">(</FONT></B>hyp<B><FONT COLOR="#663300">,</FONT></B> cov<B><FONT COLOR="#663300">,</FONT></B> lik<B><FONT COLOR="#663300">,</FONT></B> x<B><FONT COLOR="#663300">,</FONT></B> y'<B><FONT COLOR="#663300">)</FONT></B>
