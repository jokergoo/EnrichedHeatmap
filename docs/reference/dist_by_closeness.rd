<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Distance by Closeness — dist_by_closeness • EnrichedHeatmap</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Distance by Closeness — dist_by_closeness"><meta property="og:description" content="Distance by Closeness"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">
    

    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">EnrichedHeatmap</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.29.4</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li class="dropdown">
  <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
    Articles
     
    <span class="caret"></span>
  </a>
  <ul class="dropdown-menu" role="menu"><li class="dropdown-header">Vignettes</li>
    <li>
      <a href="../articles/EnrichedHeatmap_intro.html">Make Enriched Heatmaps</a>
    </li>
    <li>
      <a href="../articles/visualize_categorical_signals.html">Visualize Categorical Signals</a>
    </li>
    <li>
      <a href="../articles/row_ordering.html">Compare row ordering methods</a>
    </li>
    <li>
      <a href="../articles/roadmap.html">Visualize Comprehensive Associations in Roadmap dataset</a>
    </li>
    <li class="divider">
    <li class="dropdown-header">More applications</li>
  </ul></li>
<li>
  <a href="../reference/index.html">Reference</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/jokergoo/EnrichedHeatmap/" class="external-link">
    <span class="fab fa-github fa-lg"></span>
     
  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->

      

      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Distance by Closeness</h1>
    
    <div class="hidden name"><code>dist_by_closeness.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Distance by Closeness</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="fu">dist_by_closeness</span><span class="op">(</span><span class="va">mat</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <dl><dt>mat</dt>
<dd><p>A numeric matrix where the distance is calculated by rows.</p></dd>


</dl></div>
    <div id="details">
    <h2>Details</h2>
    <p>For two rows in the matrix, assume x_1, x_2, ..., x_n1 are the column index of none-zero values in row 1
and y_1, y_2, ... y_n2 are the column index for non-zero values in row 2, 
the distance between the two rows based on the closeness is calculated as:</p>
<p></p><div class="sourceCode"><pre><code>
    d_closeness = sum_i sum_j(|x_i - y_j|) / (n_1*n_2)  </code></pre></div>

    </div>
    <div id="value">
    <h2>Value</h2>
    

<p>A <code><a href="https://rdrr.io/r/stats/dist.html" class="external-link">dist</a></code> object.</p>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>Zuguang Gu &lt;z.gu@dkfz.de&gt;</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="va">x1</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">1</span>, <span class="fl">1</span>, <span class="fl">1</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">x2</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">1</span>, <span class="fl">1</span>, <span class="fl">1</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">x3</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">1</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">1</span>, <span class="fl">1</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span>, <span class="fl">0</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">m</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/cbind.html" class="external-link">rbind</a></span><span class="op">(</span><span class="va">x1</span>, <span class="va">x2</span>, <span class="va">x3</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/stats/dist.html" class="external-link">dist</a></span><span class="op">(</span><span class="va">m</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>          x1       x2</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> x2 1.414214         </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> x3 1.414214 1.414214</span>
<span class="r-in"><span><span class="fu">dist_by_closeness</span><span class="op">(</span><span class="va">m</span><span class="op">)</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span>          1        2</span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 2 1.222222         </span>
<span class="r-out co"><span class="r-pr">#&gt;</span> 3 2.222222 1.888889</span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by Zuguang Gu.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.0.7.</p>
</div>

      </footer></div>

  


  <style>nav[data-toggle='toc'] .nav .nav {display: block;}</style></body></html>

