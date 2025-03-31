/**
 * @file
 * Set Bazaarvoice Rating Summary drupal variables to js.
 */

(function ($, Drupal, drupalSettings) {
  "use strict";

  Drupal.behaviors.bazaarvoiceHostedRatingSummary = {
    attach: function (context, settings) {
      if (drupalSettings.bazaarvoiceRatingSummary.apiVersion === 'bv') {
        return;
      }
      $BV.configure('global', { productId : drupalSettings.bazaarvoiceRatingSummary.productid });
    }
  };

})(jQuery, Drupal, drupalSettings);
