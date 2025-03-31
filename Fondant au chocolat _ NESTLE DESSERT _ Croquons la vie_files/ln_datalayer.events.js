/**
 * @file
 * Datalayer global events.
 */

(function ($, Drupal, once, drupalSettings) {

  'use strict';

  Drupal.behaviors.ln_datalayer_api_events = {
    attach: function (context, settings) {
      if (drupalSettings.ln_datalayer && drupalSettings.ln_datalayer.events) {
        for (let event_name in drupalSettings.ln_datalayer.events) {
          let event = drupalSettings.ln_datalayer.events[event_name];
          if(event.hasOwnProperty('event')){
            if(!event.hasOwnProperty('module_name')){
              event.module_name = drupalSettings.ln_datalayer.data.module_name;
            }
            if(!event.hasOwnProperty('module_version')){
              event.module_version = drupalSettings.ln_datalayer.data.module_version;
            }
            if( event.hasOwnProperty('content_id') && event.content_id == ''){
              event.content_id = drupalSettings.ln_datalayer.data.content_id;
            }
            if( event.hasOwnProperty('content_name') && event.content_name == ''){
              event.content_name = drupalSettings.ln_datalayer.data.content_name;
            }
            if( event.hasOwnProperty('recipe_id') && event.recipe_id == ''){
              event.recipe_id = drupalSettings.ln_datalayer.data.content_id;
            }
            if( event.hasOwnProperty('recipe_name') && event.recipe_name == ''){
              event.recipe_name = drupalSettings.ln_datalayer.data.content_name;
            }
            if( event.hasOwnProperty('item_id') && event.item_id == ''){
              event.item_id = drupalSettings.ln_datalayer.data.content_id;
            }
            if( event.hasOwnProperty('item_name') && event.item_name == ''){
              event.item_name = drupalSettings.ln_datalayer.data.content_name;
            }
            window.dataLayer = window.dataLayer || [];
            window.dataLayer.push(event);
          }
        }

        delete drupalSettings.ln_datalayer.events;
      }
    }
  };

  Drupal.behaviors.ln_datalayer_global_events = {
    attach: function (context, settings) {
      //CTA link event
      $(once('ln_datalayer_cta_click', '.field--name-field-c-link a')).click(function(ev){
        let $link = $(this);

        let domain = '';
        try{
          //relative links
          let url = $link.attr('href');
          if(url.startsWith('/')){
            url = document.location.origin + url;
          }
          domain = new URL(url).hostname;
        }catch(e){}

        window.dataLayer = window.dataLayer || [];
        window.dataLayer.push({
          'event': 'cta_click',
          'event_name': 'cta_click',
          'link_classes': $link.attr('class'),
          'link_domain': domain,
          'link_id': $link.attr('id'),
          'link_text': $link.text(),
          'link_url': $link.attr('href'),
          'module_name': drupalSettings.ln_datalayer.data.module_name,
          'module_version': drupalSettings.ln_datalayer.data.module_version,
        });
      });

      //Outbound links
      $(once('ln_datalayer_outbound_link', "a[target='_blank']")).click(function(ev){
        let $link = $(this);

        let domain = '';
        try{
          //relative links
          let url = $link.attr('href');
          if(url.startsWith('/')){
            url = document.location.origin + url;
          }
          domain = new URL(url).hostname;
        }catch(e){}

        window.dataLayer = window.dataLayer || [];
        window.dataLayer.push({
          'event': 'outbound_link',
          'event_name': 'outbound_link',
          'link_classes': $link.attr('class'),
          'link_domain': domain,
          'link_id': $link.attr('id'),
          'link_url': $link.attr('href'),
          'outbound (boolean)': true,
          'module_name': drupalSettings.ln_datalayer.data.module_name,
          'module_version': drupalSettings.ln_datalayer.data.module_version,
        });
      });

      //Datalayer GA4 event - drupal_upload_submit
      $(once('ln_datalayer_input', '.form-managed-file')).click(function(ev){
        window.dataLayer = window.dataLayer || [];
        window.dataLayer.push({
          'event': 'drupal_upload_submit',
          'event_name': 'drupal_upload_submit',
          'module_name' : drupalSettings.ln_datalayer.data.module_name,
          'module_version' : drupalSettings.ln_datalayer.data.module_version,
        });
      });
      //Datalayer GA4 event - drupal_upload_error
      if ($('.form-type-managed-file').siblings(".alert-wrapper").length) {
        window.dataLayer = window.dataLayer || [];
        window.dataLayer.push({
          'event': 'drupal_upload_error',
          'event_name': 'drupal_upload_error',
          'error_name': $(".alert-wrapper .item-list").text(),
          'module_name': drupalSettings.ln_datalayer.data.module_name,
          'module_version': drupalSettings.ln_datalayer.data.module_version,
        });
      }
    }
  };
})(jQuery, Drupal, once, drupalSettings);
