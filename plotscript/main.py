def acmodel2(t, p):
    r, b, rt, bt = p
    dp = [0, 0, 0, 0]
    lambda_a = 1.0
    lambda_t = 1.0
    p_t = 0.5
    p_c = 0.5
    k = 10.0
    flag_r_0 = 1.0 if r > 0 else 0.0
    flag_b_0 = 1.0 if b > 0 else 0.0
    flag_bt_0 = 1.0 if bt > 0 else 0.0
    flag_rt_0 = 1.0 if rt > 0 else 0.0
    flag_r_rt_0 = 1.0 if r+rt > 0 else 0.0
    flag_b_bt_0 = 1.0 if b+bt > 0 else 0.0
    # R_INDEX = 0;
    dp[0] = -flag_r_0 * lambda_a * k * r / (r + bt) * p_t * r - flag_r_0 * lambda_a * k * r / (r + bt) * p_t * rt \
            + flag_r_rt_0 * lambda_t * rt - flag_r_0 * lambda_t * r / (r + rt) * rt \
            + flag_b_0 * lambda_t * (b / (b + bt)) * bt \
            + flag_bt_0 * lambda_t * (bt / (b + bt)) * bt \
            + flag_rt_0 * lambda_a * p_c * k * rt / (rt + b) * (b + bt)

    # B_INDEX = 1;
    dp[1] = -flag_b_0 * lambda_a * k * b / (b + rt) * p_t * b \
            - flag_b_0 * lambda_a * k * b / ((b + rt)) * p_t * bt \
            + flag_b_bt_0 * lambda_t * bt \
            + flag_r_0 * lambda_t * r / (r + rt) * rt \
            - flag_b_0 * lambda_t * b / (b + bt) * bt \
            + flag_rt_0 * lambda_t * (rt / (r + rt)) * rt \
            + flag_bt_0 * lambda_a * p_c * k * (bt / (r + bt)) * (r + rt)

    # RT_INEDX = 2;
    dp[2] = -flag_rt_0 * lambda_a * p_c * k * rt / (rt + b) * (b + bt) \
            - flag_r_rt_0 * lambda_t * rt \
            - flag_rt_0 * lambda_t * (rt / (r + rt)) * rt \
            + flag_r_0 * lambda_a * k * (r / (r + bt)) * p_t * r \
            + flag_r_0 * lambda_a * k * (r / (r + bt)) * p_t * rt

    # BT_INDEX = 3;
    dp[3] = -flag_bt_0 * lambda_a * k * (bt / (r + bt)) * p_c * (r + rt) \
            - flag_b_bt_0 * lambda_t * bt \
            - flag_bt_0 * lambda_t * (bt / (b + bt)) * bt \
            + flag_b_0 * lambda_a * k * b / (b + rt) * p_t * b \
            + flag_b_0 * lambda_a * k * b / (b + rt) * p_t * bt

    return dp

with open('/Users/loreti/Desktop/DATA/rb_10__b_.data','r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))


