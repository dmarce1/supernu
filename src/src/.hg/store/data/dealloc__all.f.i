        c  �       ����������$5�|cg2Eò��S�r�            x��S�n� ��+�[S������G���@d����M6Z��r����l��%�9h�С0ƵۢM�ӿ�H9�=�8��u8k�ͻ*�������f;&�PO���m���%�H���2��>B�G��ھqJy�ojt�Xg�8�hǽ�" �a��::��v�Tr'r=���ve��Ұ �X=�B[m�=�JB$�k�`a������ ��� _	�f���v�GG+��Q��*��������g��ӄ�|	.D>�r�j��>� es�e�Uw�[�/d��9�-���R�I�u��Ó�L�l���iŴy�])�'�yԹ�X�&�8�~����G?N������w    c     *  I           �����D�t���ig�M&EV�Q^��            x�c``h```�� ���*�Aiq�BzbqzQfJn~
 �2	8    �     �  �         ������`>:�uBY�z���Wڽ�            x�u�Aj�0EE�����ޔ(�:g1�4Q�J�0���t�k���Ե@3O��� ���Ve�Բ�5i�h<�h+����z~�kc��9F��i����>�a'W{���O`4M�#�};��@x>?�{�
:۱�hf{�m%c�4��C��6�4���>�`�,`t��?,����BC��5n�o���D��M��y���`t�,��c������""⩣������!hW^x�V�Ĺ���~ �R�Q    �     :  .         ������N^�Z����(o��            x�c``if c+0HIM���ON,I�HO,�OJ�OM���3�s��4�`jnf�& ���    �        .      #   ����*���j����g���^s               �   �         IMPLICIT NONE
    �             +   ����R�c�7���H#�< �w�s�*              �  �        �     x  �      -   ����Vd�>�� �ZnTЊ!�(Q            x�c``�c`````�� )��99�ɉ%����y��ɩyť�: ^QbQ��biraՒ��Z VT��[�f�%��P^����X��g�礦�hr���{�nSǣ�1-9?5M� ޔAg    d     f  �      0   ����L_o��|�N5T�Sj9��ų            x�c``����l���P�������X�^�����¥� )��99�ɉ%�@���Ĝb�<L%'hf�id�d�ڂHMlZ���JRs22�K����M. ��+4    �        �      @   ����ل��Ζjɭ��gVRq@               �   �         implicit none
    �     �        ]   ����R~��������Q�y��䴐̱            x�c``:���l���0YRRsr�KR5����st@Trb�&�XAf�FfnA��-�4�D�P��Sl�S��[��Y\椔�hra����ʀx)�y�X5�����&�4��B0�9Vř�8��bzb.؈�@<h# �T}    v     .     	   n   	����K��p����&� 9�w���            x�c``1c``Ib``PR ���Ĝ���ĒT�������]�� �z�    �     -  )   
   q   
����=�F�g�K�������G^��              ,  ,   !      !write(*,*) 'here 1.5 ...'
    �     S  |      z   �����Ej����k��j���[            x�c``�d cHIM���ON,I�HO,�O-JLII�+N���[ Ļ�j���-�Lׁ�I9pVX�uf��nKEq~iQr�& �2(�    $     +  �      �   ����^��#���{�&~9���#              *  *         deallocate(gas_curvcent)
    O     4  �      �   ����F���},NK�b{jC�Q�Ž              �  �   (      deallocate(gas_siggrey,gas_fcoef)
    �     d  �      �   ����~�������&�SC_�wz            x�c``~�������� )��99�ɉ%����ɉ:P�(�8G��*�lXejnf	����$&�&f���4���V�L���X�L-JLII�+�� .3p    �     %  �      �   ����λa�04&@�҅CY1��u              �  	         deallocate(gas_wl)
         6  �      �   ����������	a>i�2���x            x�c``�d c+0HIM���ON,I�HO,�O��,с1R+�̼��8#�B� pb    B     4  �      �   ����ɽ]���o���=��Ӆ6�;            x�c``����|���AORRsr�KR5���SSRt@���\#I� ��    v     �  3      �   ������Yz����/���9���@�
            x�c``fd``.f``X� �i����� �@S!%51''?9�$U#=�8�,1��H�*I�-��,.sR�s4��2 9��(�8G�.��T$�.��30� 
эi�\4������I`���sR�sPxE X&20�*-�B�Ǌ3�ӋR+qX������C�����L���*"RWjEq~iQr*,��s5� /T��    7     �  a      �   ����b�]A/�J���T(��:;            x�u��� �s����b)��d%Y� ��g����kfv �/����'��j�����]{U}��f^�����dI��nS��0�դt�N����n.'o������2�,	��r��/iyr��f/    �     "  w      �   ����P����;z�!�ܙw��e;��               �   �         use inputstrmod
    �     7  �      �   �����v�}3�}x�Iiƶ]�              4  O   +      if(impi==impi0) deallocate(gas_rhob)
         K  �      �   ����-{�����W��*
��C�0            x�c``Kb``����С �i����� �@S!%51''?9�$U���(�,5�(3=�D��X�9�iĪ�M,.�� �r.7    \     7  �      �   ����%U��R|�K`2�4�9����              �  �   +      deallocate(gas_eraddens,gas_luminos)
    �     N  �      �   ����`~`�!�Ύ{�?������            x�c``g``�```�T ��4��܂L[[i����������X����X_��Sl�b���dd�hr��G��0  &�2    �     9  |      �   ����(^$-�e��$
8sc�)Kh��              a  �   -      if(impi==impi0) deallocate(gas_capgam)
    	     1  �      �   �������k]�H�/�U�haf            x�c``}� �:
`��������X����X�Z�����Y���A9�\ �k    	K     ;  �     	   ����.'����w�EE���tS�y�^            x�c``Sf cC0�L���-ȴ���
)��99�ɉ%�E%�%9�E�y)�EE�\ ���    	�     k          ����$��1"7C*J�n�+J�'?�;�            x�c``�� �b
`PZ����WPZR�X���������!@��i����� �@/1/EO//�D/3/>/�(��$��DS!%51''?9�$U���$�$���2��H� ��".    	�     B       p   ������K"�>潒:�s]��6S�            x�c``��������`� �i����� �@S!%51''?9�$U#=�8�83=�(�R������ �x3    
3     E  <     �   ������K�W�Ļ���(����Zy            x�c``����r���!XRRsr�KR5���S�SRR�5��J��f���@�)�e0f^i�& `u    
x     (  X     �   �����U���}�ha��
�Q���              x  x         deallocate(gas_sigbb)
    
�     I  K      �   ������#����o�%?�S'0�            x�c``~���"���`� )��99�ɉ%����)�: FIj.��f妖d��'hr10�V ��A �
    
�     =     !  �   !����<���^����i�_�^2���            x�c``~�������`� )��99�ɉ%����)�: FIj.���Z���R\�X�� %F    &       �   "  �   "�����*��p�p&���]:�n�u�T              �  �        2     H  V   #  �   #������O���k$�����~�A˓
            x�c``Nc``^� ,�@�d�(�AJjbNN~rbI�Fzbq|jnfIAQ~���_������������  ��    z     J  �   $  �   $����m��H�(v�j|H�E%�            x�c``.e``^���P� )��99�ɉ%����y��ɩyťŚ\X�+��t@�J�
���8�:�<� B,�    �       p   %     %������>;���D��T_@Z�;'���              �  �        �     8  ?   &     &������� �!���q�L5����-                x   ,      if(impi==impi0) deallocate(str_xleft)
         X  5   '     '�������#Ր��t�y'���w�Z            x�c``�e``Nc``�Q ��4��܂L[[i����������X����X_��Sl������ħ����W���^�Z���.(J-�,�� ��&�    `     :  #   (  H   (����ݑr�կ�
�TC._r�9)              �  �   .      if(impi==impi0) deallocate(gas_siggrey)
    �       �   )  M   )�����8����CSSd�Xcg              -  Y        �       �   *  N   *����M�H�(���"�2 ���              a  �        �     U  �   +  Q   +����?��R���]���b�S�s�d��            x�c``Qd``	e �̐������SZ���¥ )��99�ɉ%�i9��9: *�L��jbW�S����_�e�����y�@= �%�           �   ,  T   ,����1�N}�M}��JBCyH�*a�f               �   �         use fluxmod
    %     �  �   -  t   -�����j�T[���7�&�Cn�g            x���M� �1�MOᲘ.�@���S$���zV#�օ�,&y�y���lv���<A)� b-<?�^34�����}3���=A��]�D�/�ϖ��B��2:oO��u�¥���y��Ȳ��v�m����� ��8�u�bW�R��?R1�G�f2��Uo���L�4s��Z=@d�    �          .  w   .������ۮ���כZҢ�;�=�               �   �         use gridmod
    �    C     0  }   /�����8}��9M���A��l
�&*            x��UM��0��+�ӆa�C/��a�C��"���l6�__�@��a} ���lK��÷%�6(� Qh����I58��#|�z�(�d�����,>~��Y?�8[���8e�E�9���x�4���4.Qe�6�@+(�h����kU� �YL�4�ÔV؀p���%Hj��Jp�`6��$L�g���� D��@h��ʞ''Gu@b��p�Ay+i�W.���B��,ԍ�<�W0�I�A���w��D��e�:�Nӽ��
�?d��.H?�_����T�z{�t:*�Ur���c.;�%�|�˲�|���?B������m�Zs���5�:��d'Y��,��U��K�¶�B�[���;�j�~2����U2JlN@s3����*�I��A�w[9ĉ{���⒐���n�����L�o�+�����_��7v��޲��n0��w�-%;�t�;3~���|�,��+�ܺ�+[XG�+o(A����*G�4��M�|����Nc����5�ϋ�y�&gM骨�M�ek�lćM-������,�L;��Z��C�FY�і�1�|��dŻu�?���3�    :     :     0  �   0������ӶC�]��B�=��Ї��              9  g   .      if(impi==impi0) deallocate(grd_capgrey)
    t     6  �   1  �   1����dJ7���~��6&z ��z�&:�            x�c``9��������� )��99�ɉ%��E)񩹙%:0Fj��W���� �V    �       �   2  Z   2������_	��<��Z2���}a��              $  I        �     +  �   3  c   3����+!�++'���Z"���              �  �         deallocate(grd_opacleak)
    �       �   4  r   4����-��H$:8 �2z-�`_D�k               �   �      j  ?        �     3     5  x   5����NI�Jg���,�K�'�>�            x�c``Nc cS0�L�H���ON,IM�H/J��LL�/��THI��#	s ��    ,       �   6  �   5�����d<"^U�Xu�/�Vy��@��              f  �      �  �        D     3  �   7  �   7   6~n�1�'Z~��J�9r�w�_&�            x�c``Nc cS0�L�H���ON,IM�H/J��LL�/��THI��#	s ��    w     �  �   8  �   8�����3��W�E�A�0
`y���            x�c``�b``h be0(-NU�-���O�Bd���D�ʖ@�0B6�(�� "͔�ۀ�; �ɉ99`��)�@f~2�xz"v��l�i9����J2�sR��A|�&�%��j�ė��@�̼2M���X��  e�N�         0  �   9     9�����N�l_�$�>-�v���Դ            x�c``z��������� )��99�ɉ%��E��9%:P:3�L� ��    2     5  �   :  x   :����Gy\ ��o�G1���R0���            x�c``�a cC0HIM���ON,I�H/*�/�Iсҙye`fqAjrbQ�& ���    g     %  z   ;  �   ;����VV���m5ti�薾��0kk              �  =         call group_dealloc
    �     U  �   <  �   9����4�	|� ��3�!
���9            x�c``����t���AWRRsr�KR5
�J��J2�sR�u@�����ļM.����@}�Ӌ
��st Tf^�& ״V    �     9  �   =  �   =����J7E�M�]'�\	A;�J��E              �  �   -      deallocate(trn_particles,trn_isvacant)
         A     >  �   9�����hΦ �"���\;��7w�H            x�c``����p��6�= CWRRsr�KR5
�J��J2�sR�u@�����ļM. �,    [     Z  �   ?  Y   ?����t��&���\huNj�;r!�            x�c``8� �
`PZ��P�_Z�����+J�K���100�B0�D>%51''?9�$U��(9>�,191�(�Ʌ��(/%���*�� �$"e    �     @     @  c   @����7s���+��7d|l�կnC            x�c``�� �b
`PZ��P���Z\�Z�������|��!*RRsr�KR5J��K��4� ~2"    �     �  �   A  �   A����*۟��0�L{۹'�IF�$�            x�=��
�0�=�)�E�Z| �P��V�c0�H���B_�g3t�3�1���q�vu'	��8�7?��=�O͡!Iܒ�W�\Kpdz;?�&��Ȋ˘�a��^��9�X�z�L״��
��.لu@9q��ᡥ��&W'�:�(���Z߇"G)��E,���>�>�    �       �   B  �   B����˔�w�nQ�8 P�-�!���              �  �   c vim: fdm=marker
    �     g  �   C  �   C�����fL�b[4����bm�LH               c   �   [*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
