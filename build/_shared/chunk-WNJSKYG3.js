import{a as P,b as S}from"/tle-finitevolume/build/_shared/chunk-22THSYFL.js";import{b as I,d as K}from"/tle-finitevolume/build/_shared/chunk-RAQ24GF6.js";function U(r){return Array.isArray(r)?r:typeof r=="number"?[y,r]:[r]}var y,c,a,h,T=I(()=>{S();y=!0,c=!1,a="skip",h=function(r,n,f,i){typeof n=="function"&&typeof f!="function"&&(i=f,f=n,n=null);let g=P(n),u=i?-1:1;l(r,void 0,[])();function l(t,C,s){let p=t&&typeof t=="object"?t:{};if(typeof p.type=="string"){let e=typeof p.tagName=="string"?p.tagName:typeof p.name=="string"?p.name:void 0;Object.defineProperty(x,"name",{value:"node ("+(t.type+(e?"<"+e+">":""))+")"})}return x;function x(){let e=[],m,o,E;if((!n||g(t,C,s[s.length-1]||null))&&(e=U(f(t,s)),e[0]===c))return e;if(t.children&&e[0]!==a)for(o=(i?t.children.length:-1)+u,E=s.concat(t);o>-1&&o<t.children.length;){if(m=l(t.children[o],o,E)(),m[0]===c)return m;o=typeof m[1]=="number"?m[1]:o+u}return e}}}});var N=I(()=>{T()});var O,b=I(()=>{N();N();O=function(r,n,f,i){typeof n=="function"&&typeof f!="function"&&(i=f,f=n,n=null),h(r,n,g,i);function g(u,l){let t=l[l.length-1];return f(u,t?t.children.indexOf(u):null,t)}}});var X={};K(X,{CONTINUE:()=>y,EXIT:()=>c,SKIP:()=>a,visit:()=>O});var d=I(()=>{b()});export{y as a,c as b,a as c,h as d,N as e,O as f,X as g,d as h};
