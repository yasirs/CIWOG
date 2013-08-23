#!/usr/bin/perl

use CGI ":all";

require 'conf.pl';

print header;

print "<html>

<body>
<form method='post' action='display.pl'>

<img border=0 src='${HTMLurl}$imgTitle'>
<br>

Web Displays<Br>
<ul>
<li>Explore gene groups. <br> Enter a Cluster ID: <input type='text' name='cid' value='6514'> <input type='submit' value='Display'>

</ul>



</ul>

Downloads<Br>
<ul>
<li><a href='exportMatrices.pl'>Download</a> cintron tabular matrices.
<li><a href='exportIntronSequences.pl'>Download</a> cintron nucleotide sequences.
</ul>


</form>
</body>


</html>";
