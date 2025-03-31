/**
 * @file
 * Shows throbber on gigya login and logout.
 */

(function ($, Drupal) {

  'use strict';

  var showThrobber = function () {
    $('body').append('<div class="ajax-progress ajax-progress-fullscreen">&nbsp;</div>');
  }

  /**
   * @type {{attach: Drupal.behaviors.clvGigyaLogin.attach}}
   */
  Drupal.behaviors.clvGigyaThrobber = {
    attach: function (context, settings) {
      gigyaHelper.addGigyaFunctionCall('accounts.addEventHandlers', {
        onLogin: showThrobber,
        onLogout: showThrobber
      });
    }
  };

})(jQuery, Drupal);
