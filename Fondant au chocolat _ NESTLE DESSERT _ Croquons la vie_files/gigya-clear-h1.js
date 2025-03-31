/**
 * @file
 * Remove all h1 that is empty.
 */

(function ($, Drupal, drupalSettings, once) {
  "use strict";

  Drupal.behaviors.clvGigyaClearH1 = {
    attach: function (context, settings) {
      if (drupalSettings.gigya.enableRaaS) {
        if (drupalSettings.gigya.raas.register) {
          drupalSettings.gigya.raas.register.onAfterScreenLoad = this.clearGigyaCaption;
        }
      }
    },
    clearGigyaCaption: function (e) {
      var containerID = e.sourceContainerID;
      if (typeof containerID !== undefined) {
        // Remove empty screenset caption.
        $(once('cleanH1', "#" + containerID + " h1.gigya-screen-caption"))
          .each(function () {
            const h1 = $(this);
            if (h1.html().replace(/\s|&nbsp;/g, "").length === 0) {
              h1.remove();
            }
          });
      }
    },
  };
})(jQuery, Drupal, drupalSettings, once);
