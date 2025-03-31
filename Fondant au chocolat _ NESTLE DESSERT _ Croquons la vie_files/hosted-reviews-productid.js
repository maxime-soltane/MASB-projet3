(function ($) {
  "use strict";

  Drupal.behaviors.bazaarvoiceHostedReviews = {
    attach: function (context, settings) {
      if(drupalSettings.bazaarvoiceReviews.apiVersion === 'bv') {
        return;
      }
      $BV.configure('global', { productId : drupalSettings.bazaarvoiceReviews.productid });
    }
  };

})(jQuery);


